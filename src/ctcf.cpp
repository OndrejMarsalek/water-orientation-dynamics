///=============================================================================
///
/// Conditional time correlation functions
///
/// author: Ondrej Marsalek
///         ondrej.marsalek@gmail.com
///
///=============================================================================


//------------------------------------------------------------------------------
// includes and namespace stuff
//------------------------------------------------------------------------------

// system includes
#include <deque>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>

// This cannot be put in gmx::, as it pulls in the new C++ headers
// that already have the namespace.
// However, `gmx_run_cmain` is the only thing in this header outside of gmx::.
#include "gromacs/commandline/cmdlineinit.h"

// keep GMX stuff isolated
namespace gmx {
#include "gromacs/legacyheaders/macros.h"
#include "gromacs/topology/topology.h"
#include "gromacs/topology/index.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/math/vec.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/fileio/trxio.h"
}

// pull in std:: members explicitly
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::deque;
using std::vector;
using std::ofstream;
using std::setw;
using std::setprecision;
using std::ios;
using std::min;

// pull in GMX stuff explicitly
using gmx::etREAL;
using gmx::etINT;
using gmx::etENUM;
using gmx::real;
using gmx::rvec;
using gmx::atom_id;
using gmx::t_topology;
using gmx::gmx_fatal;
using gmx::read_top;
using gmx::matrix;
using gmx::rvec_add;
using gmx::rvec_sub;
using gmx::copy_rvec;
using gmx::svmul;
using gmx::invnorm;
using gmx::norm;
using gmx::cos_angle;
using gmx::read_first_x;
using gmx::read_next_x;
using gmx::close_trj;
using gmx::t_pargs;
using gmx::t_filenm;
using gmx::efTRX;
using gmx::efTPR;
using gmx::efNDX;
using gmx::efDAT;
using gmx::t_pbc;
using gmx::gmx_rmpbc_t;
using gmx::gmx_rmpbc_init;
using gmx::gmx_rmpbc;
using gmx::output_env_t;
using gmx::t_trxstatus;
using gmx::gmx_ffopen;
using gmx::gmx_ffclose;


// list of possible types of distance coordinate
enum ecoord {
    coordR,
    coordZ
};


// order of Legendre polynomial
enum ePorder {Porder1, Porder2};


//------------------------------------------------------------------------------
// the 2nd Legendre polynomial
//------------------------------------------------------------------------------
template <class T>
T P2(T x) {
    return 1.5*x*x - 0.5;
}


//-------------------------------------------------------------------------------
// rvec cache class
//-------------------------------------------------------------------------------
class Cache {

private:
    deque<rvec*> data;
    unsigned int num;

public:
    Cache(const unsigned int num);
    void push_back(const rvec* x);
    rvec* peek(const unsigned int i) const;
    void pop_front();
    int size() const;
};

Cache::Cache(const unsigned int num):
    num(num)
{}

void Cache::push_back(const rvec* x) {
    // x must point to num rvec
    rvec* n = NULL;
    n = new rvec[num];
    data.push_back(n);
    for (unsigned int i=0; i<num; ++i) { // copy current positions to cache
        copy_rvec(x[i], n[i]);
    }
}

void Cache::pop_front() {
    delete[] data.front();
    data.pop_front();
}

rvec* Cache::peek(const unsigned int i) const {
    return data[i];
}

int Cache::size() const {
    return data.size();
}


//------------------------------------------------------------------------------
//  conditional timeseries
//------------------------------------------------------------------------------
template <class T>
class CTS {

private:
    vector<T> data;
    double r_max;
    unsigned int n_r;
    double dr;
    unsigned int n_tau;
    ePorder Porder;

public:
    CTS(const double r_max, const unsigned int n_r, const unsigned int n_tau,
        const ePorder Porder);
    void write(const string filename, const vector<double> &taus) const;
    inline void set(const unsigned int i, const unsigned int j, const T value);
    T get(const unsigned int i, const unsigned int j) const;
    void add(const double d_ref, const Cache& c);
};

template <class T>
CTS<T>::CTS(const double dr, const unsigned int n_r, const unsigned int n_tau,
            const ePorder Porder):
    r_max(dr * n_r),
    n_r(n_r),
    dr(dr),
    n_tau(n_tau),
    Porder(Porder)
{

    // allocate and zero storage
    data.resize(n_r * n_tau, 0);

}

template <class T>
void CTS<T>::write(const string filename, const vector<double> &taus) const {

    // open file with GMX wrapper to get backups
    FILE* f_out = gmx_ffopen(filename.c_str(), "w");

    // write Legendre polynomial order
    if (Porder == Porder1) {
        fprintf(f_out, "%12.6f", 1.0);
    } else if (Porder == Porder2) {
        fprintf(f_out, "%12.6f", 2.0);
    }

    // write delay time tau
    for (unsigned int i=0; i<n_tau; ++i) {
        fprintf(f_out, "%12.6f ", taus[i] - taus[0]);
    }
    fprintf(f_out, "\n");

    int c = 0;
    for (unsigned int j=0; j<n_r; ++j) {

        // write r bin center
        fprintf(f_out, "%12.6f ", (j+0.5) * dr);

        // write data for this r bin
        for (unsigned int i=0; i<n_tau; ++i) {
            fprintf(f_out, "%12.6f ", data[c]);
            c++;
        }
        fprintf(f_out, "\n");
    }

    gmx_ffclose(f_out);
}

template <class T>
void CTS<T>::add(const double d_ref, const Cache& c) {

    if ((0 < d_ref) && (d_ref < r_max)) {
        const unsigned int bin_r = int(d_ref / dr);
        const unsigned int i_base = bin_r * n_tau;
        const T* vec_base = c.peek(0)[0];
        if (Porder == Porder1) {
            for (unsigned int i=0; i<n_tau; ++i) {
                data[i_base + i] += cos_angle(vec_base, c.peek(i)[0]);
            }
        } else if (Porder == Porder2) {
            for (unsigned int i=0; i<n_tau; ++i) {
                data[i_base + i] += P2(cos_angle(vec_base, c.peek(i)[0]));
            }
        }
    }

}


//------------------------------------------------------------------------------
// the worker function, process the whole trajectory
//------------------------------------------------------------------------------
void process_trj(const char* fn_trj, const char* fn_ndx, const char* fn_top,
                 const char* fn_oh, const char* fn_dm,
                 const double bw_r, double r_max, double n_tau, const int stride,
                 const ePorder Porder, const ecoord coord,
                 const output_env_t output_env) {

    real t0, t1;           // time
    rvec* x;               // particle positions of current frame
    rvec* x_ref;           // particle positions from cache for reference
    rvec vec_OH1;          // 1st OH bond vector
    rvec vec_OH2;          // 2nd OH bond vector
    rvec vec_DM;           // dipole moment vector
    rvec dx_O;             // relative position
    double d_ref = 0.0;    // reference distance; suppress uninitialized warning
    vector<double> taus;   // values of tau for labelling

    // read the topology file
    t_pbc pbc;                 // PBC
    int ePBC;                  // type of PBC
    t_topology* top = read_top(fn_top, &ePBC);

    // read one group from the index file - water molecules, order O H H
    atom_id* index_w;   // index of water atoms
    int index_size_w;   // size of water index group
    char* grpname_w;    // name of water index group
    cerr << "\nSelect the water group (assumed ordering: O H H).\n";
    get_index(&top->atoms, fn_ndx, 1, &index_size_w, &index_w, &grpname_w);
    if ((index_size_w % 3) != 0) {
        gmx_fatal(FARGS, "Number of index elements (%i) not a multiple of 3, "
                         "these can not be water molecules.\n", index_size_w);
    }
    const unsigned int num_w = index_size_w / 3;   // number of water molecules

    // read one group from the index file - solute
    atom_id* index_s;   // index of solute atoms
    int index_size_s;   // size of solute index group
    char* grpname_s;    // name of solute index group
    cerr << "\nSelect the solute group (one atom).\n";
    get_index(&top->atoms, fn_ndx, 1, &index_size_s, &index_s, &grpname_s);
    cerr << "\n";

    // only one atom allowed at the moment
    if (index_size_s != 1)
        gmx_fatal(FARGS, "Only one solute atom allowed at the moment.\n");

    // store the solute atom index
    const int i_ref = index_s[0];         // solute atom index

    // read first frame and keep trajectory status
    t_trxstatus* trxstatus;
    matrix box;
    const int num_atoms = read_first_x(output_env, &trxstatus, fn_trj, &t0, &x, box);
    if (!num_atoms) {
        gmx_fatal(FARGS, "Could not read coordinates from file '%s'.\n", fn_trj);
    }

    // Default r_max is such that ALL water molecules fall under r_max.
    // That means the sum of contributions at all distances will match the unresolved total ACF.
    double r_max_box = -1.0;    // silence warnings
    if (coord == coordR) {
        // half the box diagonal
        r_max_box = 0.5 * sqrt(pow(box[0][0], 2) + pow(box[1][1], 2) + pow(box[2][2], 2));
    } else if (coord == coordZ) {
        // the whole box in the z direction
        r_max_box = box[2][2];
    }
    // Set r_max either to cover the whole box times a small fudge factor for possible pressure coupling,
    // or the value requested by the user.
    if (r_max < 0.0) {
        r_max = 1.05 * r_max_box;
    } else {
        if (r_max < r_max_box) {
            printf("\nWARNING: r_max smaller than half the box diagonal - analysis will likely miss some water molecules.");
        }
    }

    // number of bins in r
    unsigned int n_r = (unsigned int)(r_max / bw_r);

    // prepare CTCF collectors
    CTS<real> cacf_O_OH(bw_r, n_r, n_tau, Porder);  // OH bonds
    CTS<real> cacf_O_DM(bw_r, n_r, n_tau, Porder);  // dipole moment vector

    // prepare cache
    Cache cache_x(num_atoms);
    Cache sample(1);
    vector<Cache> cache_OH1(num_w, sample);
    vector<Cache> cache_OH2(num_w, sample);
    vector<Cache> cache_DM(num_w, sample);

    // initialize PBC removal
    gmx_rmpbc_t gpbc = gmx_rmpbc_init(&top->idef, ePBC, num_atoms);

    // prepare time labels array
    taus.resize(n_tau);

    // loop over frames
    unsigned int frame = 0;
    t1 = t0;
    do {

        // Make molecules whole, as `mdrun` does not guarantee whole molecules.
        gmx_rmpbc(gpbc, num_atoms, box, x);

        // update PBC with current box dimensions
        set_pbc(&pbc, ePBC, box);

        // store current positions in cache
        // (will be needed for R)
        cache_x.push_back(x);

        // trim position cache if needed
        if (cache_x.size() > n_tau) {
            cache_x.pop_front();
        }
        // loop over all water molecules
        for (unsigned int i=0; i<num_w; i++) {

            // get the OH bond and dipole moment vectors ...
            // (can use just rvec_sub here, as molecules are whole)
            // ... normalize them ...
            const int ai = index_w[3*i];     // O atom
            const int aj = index_w[3*i+1];   // H1 atom
            const int ak = index_w[3*i+2];   // H2 atom
            rvec_sub(x[ai], x[aj], vec_OH1);             // vec_OH1 is the 1st OH bond vector
            svmul(invnorm(vec_OH1), vec_OH1, vec_OH1);
            rvec_sub(x[ai], x[ak], vec_OH2);             // vec_OH2 is the 2nd OH bond vector
            svmul(invnorm(vec_OH2), vec_OH2, vec_OH2);
            rvec_add(vec_OH1, vec_OH2, vec_DM);          // vec_DM is in the direction of the dipole moment, but different length
            svmul(invnorm(vec_DM), vec_DM, vec_DM);

            // ... and cache them
            cache_OH1[i].push_back(&vec_OH1);
            cache_OH2[i].push_back(&vec_OH2);
            cache_DM[i].push_back(&vec_DM);

            // trim cache if needed
            if (cache_OH1[i].size() > n_tau) {    // all the sizes should be the same
                cache_OH1[i].pop_front();
                cache_OH2[i].pop_front();
                cache_DM[i].pop_front();
            }

            // starts of individual contributing segments are `stride` apart
            if ((frame > n_tau-2) && (frame % stride) == 0) {

                // Calculate the reference distance based on positions at the start of this segment.
                x_ref = cache_x.peek(0);
                if (coord == coordR) {
                    pbc_dx(&pbc, x_ref[i_ref], x_ref[ai], dx_O);
                    d_ref = norm(dx_O);
                } else if (coord == coordZ) {
                    d_ref = x_ref[ai][2];
                }

                // orientation CACF, using the cache of OH and DM vectors
                cacf_O_OH.add(d_ref, cache_OH1[i]);
                cacf_O_OH.add(d_ref, cache_OH2[i]);
                cacf_O_DM.add(d_ref, cache_DM[i]);
            }

        }

        // store time labels
        if (frame < n_tau) {
            taus[frame] = t1;
        }

        frame++;

    } while (read_next_x(output_env, trxstatus, &t1, x, box));
    close_trj(trxstatus);
    cerr << "\nDone with trajectory.\n";

    // write out data, pass in time labels for plotting
    cacf_O_OH.write(fn_oh, taus);
    cacf_O_DM.write(fn_dm, taus);

}


//------------------------------------------------------------------------------
// process all user input and prepare for analysis
//------------------------------------------------------------------------------
int run(int argc, char *argv[]) {

    output_env_t output_env;
    ePorder Porder = Porder1;
    static const char *Porder_name[] = {NULL, "1", "2", NULL};
    ecoord coord = coordR;
    static const char *coord_name[] = {NULL, "r", "z", NULL};

    // program description
    // TODO: write
    const char *desc[] = {
        "There will be a description here.",
    };
    const int n_desc = asize(desc);

    // default values of parameters
    real bw_r = 0.05;
    real r_max = -1.0;    // will be set based on box if not set by user
    int n_tau = 1000;
    int stride = 50;

    // specify command line options
    t_pargs pa[] = {
        {"-bin-r",       FALSE, etREAL, {&bw_r},       "Binwidth in r (nm)"},
        {"-r-max",       FALSE, etREAL, {&r_max},      "Maximum value of r (nm), <0 means automatic"},
        {"-n-tau",       FALSE, etINT,  {&n_tau},      "Number of time points of time correlation functions"},
        {"-tau-restart", FALSE, etINT,  {&stride},     "Stride of independent starts of correlation functions"},
        {"-P",           FALSE, etENUM, {Porder_name}, "Legendre polynomial order, defaults to 1"},
        {"-coord",       FALSE, etENUM, {coord_name},  "Coordinate r to use as condition"},
    };
    const int n_pargs = asize(pa);

    t_filenm fnm[] = {
        {efTRX, "-f", NULL, ffREAD},
        {efTPR, NULL, NULL, ffREAD},
        {efNDX, NULL, NULL, ffOPTRD},
        {efDAT, "-dm", "ctcf-DM", ffOPTWR},
        {efDAT, "-oh", "ctcf-OH", ffOPTWR},
    };
    const int n_file = asize(fnm);

    // parse common arguments based on the above
    if (!parse_common_args(&argc, argv, PCA_CAN_TIME,
                           n_file, fnm,
                           n_pargs, pa,
                           n_desc, desc,
                           0, NULL,
                           &output_env)) {
        return 0;
    }

    // set the Legendre polynomial order
    if (Porder_name[0][0] == '1') {
        Porder = Porder1;
    } else if (Porder_name[0][0] == '2') {
        Porder = Porder2;
    } else {
        gmx_fatal(FARGS, "Unexpected order of Legendre polynomial");
    }

    // set type of coordinate based on command line input
    if (coord_name[0][0] == 'r') {
        coord = coordR;
    } else if (coord_name[0][0] == 'z') {
        coord = coordZ;
    } else {
        gmx_fatal(FARGS, "Unexpected type of distribution function.");
    }

    // get filenames
    const char* fn_trj = ftp2fn(efTRX, n_file, fnm);
    const char* fn_top = ftp2fn(efTPR, n_file, fnm);
    const char* fn_ndx = ftp2fn_null(efNDX, n_file, fnm);
    const char* fn_dm  = opt2fn("-dm", n_file, fnm);
    const char* fn_oh  = opt2fn("-oh", n_file, fnm);

    // process the trajectory
    process_trj(fn_trj, fn_ndx, fn_top,
                fn_oh, fn_dm,
                bw_r, r_max, n_tau, stride,
                Porder, coord,
                output_env);

    return 0;
}


//------------------------------------------------------------------------------
// program entry point
//------------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    return gmx_run_cmain(argc, argv, &run);
}

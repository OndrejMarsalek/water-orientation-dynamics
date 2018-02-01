///=============================================================================
///
/// Orientation-resolved distribution functions
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

// this cannot be put in gmx::, as it pulls in the new C++ headers
// that already have the namespace.
// However, `gmx_run_cmain` is the only thing in this header outside of gmx::.
#include "gromacs/commandline/cmdlineinit.h"

// keep legacy GMX stuff isolated
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
using gmx::gmx_angle;
using gmx::read_first_x;
using gmx::read_next_x;
using gmx::close_trj;
using gmx::t_pargs;
using gmx::t_filenm;
using gmx::efTRX;
using gmx::efTPR;
using gmx::efNDX;
using gmx::efDAT;
using gmx::output_env_t;
using gmx::t_trxstatus;
using gmx::t_pbc;
using gmx::gmx_rmpbc_t;
using gmx::gmx_rmpbc_init;
using gmx::gmx_rmpbc;
using gmx::det;
using gmx::gmx_ffopen;
using gmx::gmx_ffclose;


// list of possible types of distance coordinate
enum ecoord {
    coordR,
    coordZ
};


//------------------------------------------------------------------------------
// orientation-resolved distribution function
//------------------------------------------------------------------------------
template <class T>
class OR_DF {

private:
    vector<T> data;
    double r_max;
    unsigned int n_r;
    double dr;
    double th_max;
    unsigned int n_th;
    double dth;
    ecoord coord;

public:
    OR_DF(const double dr, const unsigned int n_r,
          const unsigned int n_th, const ecoord coord);
    void write(const string filename) const;
    void add(const double d_ref, const double theta, const double irt);
    void normalize(const double n_factor);
};

//------------------------------------------------------------------------------
// construct and initialize the OR-DF collector
//------------------------------------------------------------------------------
template <class T>
OR_DF<T>::OR_DF(const double dr, const unsigned int n_r,
                const unsigned int n_th, const ecoord coord):
    r_max(dr * n_r),
    n_r(n_r),
    dr(dr),
    th_max(M_PI),
    n_th(n_th),
    dth(th_max / n_th),
    coord(coord) {

    data.resize(n_r * n_th, 0);

}

//------------------------------------------------------------------------------
// add a contribution to the OR-DF
//------------------------------------------------------------------------------
template <class T>
void OR_DF<T>::add(const double d_ref, const double theta, const double inc) {

    if ((0 < d_ref) && (d_ref < r_max) && (theta >= 0.0) && (theta <= M_PI)) {

        // determine bin indices
        const unsigned int bin_r = int(d_ref / dr);
        const unsigned int bin_th = int(theta / dth);

        // add increment
        data[bin_r + bin_th*n_r] += inc;
    }

}

//------------------------------------------------------------------------------
// normalize the OR-DF when all data is collected
//------------------------------------------------------------------------------
template <class T>
void OR_DF<T>::normalize(const double n_factor) {

    vector<T> normr;
    vector<T> normth;

    normr.resize(n_r, 0);
    normth.resize(n_th, 0);

    // prepare volume element factors in r
    if (coord == coordR) {
        real r = 0.5 * dr;
        for (unsigned int j=0; j<n_r; ++j) {
            normr[j] = (4.0/3.0) * M_PI * dr * (3*r*r + 0.25*dr*dr);
            r += dr;
        }
    } else if (coord == coordZ) {
        real r = 0.5 * dr;
        for (unsigned int j=0; j<n_r; ++j) {
            normr[j] = dr;
            r += dr;
        }
    }

    // prepare volume element factors in theta
    real th = 0.5 * dth;
    for (unsigned int i=0; i<n_th; ++i) {
        normth[i] = dth * 0.5 * sin(th);
        th += dth;
    }

    // apply normalization
    for (unsigned int i=0; i<n_th; ++i) {
        for (unsigned int j=0; j<n_r; ++j) {
            data[i*n_r + j] /= normth[i] * normr[j] * n_factor;
        }
    }

}

//------------------------------------------------------------------------------
// write the OR-DF to a file
//------------------------------------------------------------------------------
template <class T>
void OR_DF<T>::write(const string filename) const {

    // open file with GMX wrapper to get backups
    FILE* f_out = gmx_ffopen(filename.c_str(), "w");

    // use the corner element to mark the type of analysis that was done
    if (coord == coordR) {
        fprintf(f_out, "%12.6f", 0.0);
    } else if (coord == coordZ) {
        fprintf(f_out, "%12.6f", 1.0);
    }

    // write r bin centers
    for (unsigned int i=0; i<n_r; ++i) {
        fprintf(f_out, "%12.6f ", (i+0.5) * dr);
    }
    fprintf(f_out, "\n");

    int c = 0;
    for (unsigned int i=0; i<n_th; ++i) {

        // write theta bin center in degrees
        fprintf(f_out, "%12.6f ", (i+0.5) * 180.0 * dth / M_PI);

        // write data
        for (unsigned int j=0; j<n_r; ++j) {
            fprintf(f_out, "%12.6f ", data[c]);
            c++;
        }
        fprintf(f_out, "\n");
    }

    gmx_ffclose(f_out);

}


//------------------------------------------------------------------------------
// the worker function, process the whole trajectory
//------------------------------------------------------------------------------
void process_trj(const char* fn_trj, const char* fn_ndx, const char* fn_top,
                 const char* fn_oh, const char* fn_dm,
                 const double bw_r, double r_max, const double bw_th,
                 const ecoord coord, const output_env_t output_env) {

    rvec* x;                // "array" of particle positions
    rvec vec_OH1;           // 1st OH bond vector
    rvec vec_OH2;           // 2nd OH bond vector
    rvec vec_DM;            // dipole moment vector
    rvec dx_ref;            // reference vector - radial vector solute-O or z coordinate
    double d_ref = 0.0;     // length of reference, init to silence warning
    real t;                 // time
    double V = 0.0;         // box volume, init to silence warning
    t_pbc pbc;              // PBC
    int ePBC;               // type of PBC
    atom_id* index_w;       // index of water atoms
    int index_size_w;       // size of water index group
    char* grpname_w;        // name of water index group
    atom_id* index_s;       // index of solute atoms
    int index_size_s = 1;   // size of solute index group
    char* grpname_s;        // name of solute index group
    double norm = 1.0;      // normalization factor, init to silence warning

    // read the topology file
    t_topology* top = read_top(fn_top, &ePBC);

    // read one group from the index file - water molecules, order O H H
    cerr << "\nSelect the water group (assumed ordering: O H H).\n";
    get_index(&top->atoms, fn_ndx, 1, &index_size_w, &index_w, &grpname_w);
    if ((index_size_w % 3) != 0) {
        gmx_fatal(FARGS, "Number of index elements (%i) not a multiple of 3, "
                         "these can not be water molecules.\n", index_size_w);
    }
    const unsigned int num_w = index_size_w / 3;   // number of water molecules

    if (coord == coordR) {
        // read one group from the index file - solute
        cerr << "\nSelect the solute group.\n";
        get_index(&top->atoms, fn_ndx, 1, &index_size_s, &index_s, &grpname_s);
        cerr << "\n";
    } else if (coord == coordZ) {
        // no solutes, just the z coordinate
        index_size_s = 1;
    }

    // read first frame and keep trajectory status
    t_trxstatus* trxstatus;
    matrix box;
    const int num_atoms = read_first_x(output_env, &trxstatus, fn_trj, &t, &x, box);
    if (!num_atoms) {
        gmx_fatal(FARGS, "Could not read coordinates from file '%s'.\n", fn_trj);
    }

    if (coord == coordR) {
        // keep r_max reasonable if too large or unset
        double r_max_box = 0.49 * min(box[0][0], min(box[1][1], box[2][2]));
        if (r_max > r_max_box) {
            printf("\nr_max larger than half of the box, setting r_max = %6.4f\n", r_max_box);
            r_max = r_max_box;
        } else if (r_max < 0.0) {
            r_max = r_max_box;
        }
    } else if (coord == coordZ) {
        // if r_max is unset or larger than z box size, set it to z box size
        if ((r_max < 0) || (r_max > box[2][2])) {
            r_max = box[2][2];
        }
    }

    // number of bins in r
    unsigned int n_r = (unsigned int)(r_max / bw_r);
    // This is the last time the explicit r_max is used. If it does not fit
    // exactly (likely), r_max in the histogram and in plotting is actually
    // determined by n_r*bw_r.

    // number of bins in theta
    unsigned int n_th = (unsigned int)(180.0 / bw_th);

    // collectors of orientation-resolved distribution functions
    OR_DF<real> or_df_DM(bw_r, n_r, n_th, coord);
    OR_DF<real> or_df_OH(bw_r, n_r, n_th, coord);

    // initialize PBC removal
    gmx_rmpbc_t gpbc = gmx_rmpbc_init(&top->idef, ePBC, num_atoms);

    unsigned long long frame = 0;

    // loop over frames
    do {

        // Make molecules whole, as `mdrun` does not guarantee whole molecules.
        gmx_rmpbc(gpbc, num_atoms, box, x);

        // update PBC with current box dimensions
        set_pbc(&pbc, ePBC, box);

        if (coord == coordR) {
            // volume for this frame
            V = det(box);
        } else if (coord == coordZ) {
            V = 1.0 / (box[0][0] * box[1][1]);
        }

        // loop over all reference atoms
        for (int j=0; j<index_size_s; ++j) {

            // loop over all water molecules
            for (unsigned int i=0; i<num_w; i++) {

                // get the OH bond and dipole moment vectors...
                // (can use just rvec_sub here, as molecules are whole)
                const int ai = index_w[3*i];          // O atom
                const int aj = index_w[3*i+1];        // H1 atom
                const int ak = index_w[3*i+2];        // H2 atom
                rvec_sub(x[aj], x[ai], vec_OH1);      // vec_OH1 is the 1st OH bond vector
                rvec_sub(x[ak], x[ai], vec_OH2);      // vec_OH1 is the 2nd OH bond vector
                rvec_add(vec_OH1, vec_OH2, vec_DM);   // vec_DM is in the direction of the dipole moment, but different length

                // get the reference vector
                if (coord == coordR) {
                    pbc_dx(&pbc, x[index_s[j]], x[ai], dx_ref);
                    d_ref = gmx::norm(dx_ref);
                } else if (coord == coordZ) {
                    copy_rvec(x[ai], dx_ref);
                    dx_ref[0] = 0.0;
                    dx_ref[1] = 0.0;
                    d_ref = dx_ref[2];
                }

                // orientation-resolved DF
                or_df_DM.add(d_ref, gmx_angle(vec_DM, dx_ref), V);
                or_df_OH.add(d_ref, gmx_angle(vec_OH1, dx_ref), V);
                or_df_OH.add(d_ref, gmx_angle(vec_OH2, dx_ref), V);
            }
        }

        frame++;

    } while (read_next_x(output_env, trxstatus, &t, x, box));
    close_trj(trxstatus);
    cerr << "\nDone with trajectory.\n";

    // normalize OR-DF by the number of frames (and possibly pairs)
    if (coord == coordR) {
        // number of frames times contributing pairs
        norm = static_cast<double>(frame) * index_size_s * num_w;
    } else if (coord == coordZ) {
        // only number of frames
        norm = frame;
    }
    or_df_DM.normalize(norm);
    or_df_OH.normalize(2 * norm);   // factor of 2 for 2 OH bonds per O atom

    // write out data
    or_df_DM.write(fn_dm);
    or_df_OH.write(fn_oh);

}


//------------------------------------------------------------------------------
// process all user input and prepare for analysis
//------------------------------------------------------------------------------
int run(int argc, char *argv[]) {

    output_env_t output_env;
    ecoord coord = coordR;
    static const char *coord_name[] = {NULL, "radial", "z", NULL};

    // program description
    // TODO: write
    const char *desc[] = {
        "There will be a description here.",
    };
    const int n_desc = asize(desc);

    // default values of parameters
    real bw_th = 2;
    real bw_r = 0.01;
    real r_max = -1.0;  // will be set based on box if not set by user

    // specify command line options
    t_pargs pa[] = {
        {"-bin-r",     FALSE, etREAL, {&bw_r},      "Binwidth in r (nm)"},
        {"-bin-theta", FALSE, etREAL, {&bw_th},     "Binwidth in theta (deg)"},
        {"-r-max",     FALSE, etREAL, {&r_max},     "Maximum value of r (nm), <0 means automatic"},
        {"-coord",     FALSE, etENUM, {coord_name}, "Distance coordinate r to use"},
    };
    const int n_pargs = asize(pa);

    t_filenm fnm[] = {
        {efTRX, "-f", NULL, ffREAD},
        {efTPR, NULL, NULL, ffREAD},
        {efNDX, NULL, NULL, ffOPTRD},
        {efDAT, "-dm", "or-df-DM", ffOPTWR},
        {efDAT, "-oh", "or-df-OH", ffOPTWR},
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

    // get filenames
    const char* fn_trj = ftp2fn(efTRX, n_file, fnm);
    const char* fn_top = ftp2fn(efTPR, n_file, fnm);
    const char* fn_ndx = ftp2fn_null(efNDX, n_file, fnm);
    const char* fn_dm  = opt2fn("-dm", n_file, fnm);
    const char* fn_oh  = opt2fn("-oh", n_file, fnm);

    // set type of coordinate based on command line input
    if (coord_name[0][0] == 'r') {
        coord = coordR;
    } else if (coord_name[0][0] == 'z') {
        coord = coordZ;
    } else {
        gmx_fatal(FARGS, "Unexpected type of distribution function.");
    }

    // process the trajectory
    process_trj(fn_trj, fn_ndx, fn_top,
                fn_oh, fn_dm,
                bw_r, r_max, bw_th,
                coord, output_env);

    return 0;
}


//------------------------------------------------------------------------------
// program entry point
//------------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    return gmx_run_cmain(argc, argv, &run);
}

// fast5_interface needs cleaning
#define BANANA 1
#include <assert.h>
#include <err.h>
#include <math.h>
#include "fast5_interface.h"
#include "scrappie_stdlib.h"
//#include "util.h"

struct _gop_data {
    const char *prefix;
    int latest;
};

typedef struct {
    //  Information for scaling raw data from ADC values to pA
    float digitisation;
    float offset;
    float range;
    float sample_rate;
} fast5_raw_scaling;

float read_float_attribute(hid_t group, const char *attribute) {
    float val = NAN;
    if (group < 0) {
        warnx("Invalid group passed to %s:%d.", __FILE__, __LINE__);
        return val;
    }

    hid_t attr = H5Aopen(group, attribute, H5P_DEFAULT);
    if (attr < 0) {
        warnx("Failed to open attribute '%s' for reading.", attribute);
        return val;
    }

    H5Aread(attr, H5T_NATIVE_FLOAT, &val);
    H5Aclose(attr);

    return val;
}

fast5_raw_scaling get_raw_scaling(hid_t hdf5file) {
    // Add 1e-5 to sensible sample rate as a sentinel value
    fast5_raw_scaling scaling = { NAN, NAN, NAN, NAN };
    const char *scaling_path = "/UniqueGlobalKey/channel_id";

    hid_t scaling_group = H5Gopen(hdf5file, scaling_path, H5P_DEFAULT);
    if (scaling_group < 0) {
        warnx("Failed to group %s.", scaling_path);
        return scaling;
    }

    scaling.digitisation = read_float_attribute(scaling_group, "digitisation");
    scaling.offset = read_float_attribute(scaling_group, "offset");
    scaling.range = read_float_attribute(scaling_group, "range");
    scaling.sample_rate = read_float_attribute(scaling_group, "sampling_rate");

    H5Gclose(scaling_group);

    return scaling;
}

raw_table read_raw(const char *filename, bool scale_to_pA) {
    assert(NULL != filename);
    raw_table rawtbl = { 0, 0, 0, NULL };

    hid_t hdf5file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (hdf5file < 0) {
        warnx("Failed to open %s for reading.", filename);
        return rawtbl;
    }
    H5Eset_auto2(H5E_DEFAULT, NULL, NULL);

    const char *root = "/Raw/Reads/";
    const int rootstr_len = strlen(root);
    ssize_t size =
        H5Lget_name_by_idx(hdf5file, root, H5_INDEX_NAME, H5_ITER_INC, 0, NULL,
                           0, H5P_DEFAULT);
    if (size < 0) {
        warnx("Failed find read name under %s.", root);
        goto cleanup1;
    }
    char *name = calloc(1 + size, sizeof(char));
    H5Lget_name_by_idx(hdf5file, root, H5_INDEX_NAME, H5_ITER_INC, 0, name,
                       1 + size, H5P_DEFAULT);

    // Create group name
    char *signal_path = calloc(rootstr_len + size + 8, sizeof(char));
    (void)snprintf(signal_path, rootstr_len + size + 8, "%s%s/Signal", root,
                   name);
    free(name);

    hid_t dset = H5Dopen(hdf5file, signal_path, H5P_DEFAULT);
    if (dset < 0) {
        warnx("Failed to open dataset '%s' to read raw signal from.",
              signal_path);
        goto cleanup2;
    }

    hid_t space = H5Dget_space(dset);
    if (space < 0) {
        warnx("Failed to create copy of dataspace for raw signal %s.",
              signal_path);
        goto cleanup3;
    }
    hsize_t nsample;
    H5Sget_simple_extent_dims(space, &nsample, NULL);
    float *rawptr = calloc(nsample, sizeof(float));
    herr_t status =
        H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rawptr);
    if (status < 0) {
        free(rawptr);
        warnx("Failed to read raw data from dataset %s.", signal_path);
        goto cleanup4;
    }
    rawtbl = (raw_table) {
    nsample, 0, nsample, rawptr};

    if (scale_to_pA) {
        const fast5_raw_scaling scaling = get_raw_scaling(hdf5file);
        const float raw_unit = scaling.range / scaling.digitisation;
        for (size_t i = 0; i < nsample; i++) {
            rawptr[i] = (rawptr[i] + scaling.offset) * raw_unit;
        }
    }

 cleanup4:
    H5Sclose(space);
 cleanup3:
    H5Dclose(dset);
 cleanup2:
    free(signal_path);
 cleanup1:
    H5Fclose(hdf5file);

    return rawtbl;
}

void write_annotated_events(hid_t hdf5file, const char *readname,
                            const event_table et, hsize_t chunk_size,
                            int compression_level) {
    assert(compression_level >= 0 && compression_level <= 9);

    // Memory representation
    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(event_t));
    if (memtype < 0) {
        warnx("Failed to create memroy representation for event table %s:%d.",
              __FILE__, __LINE__);
        goto clean1;
    }
    H5Tinsert(memtype, "start", HOFFSET(event_t, start), H5T_NATIVE_UINT64);
    H5Tinsert(memtype, "length", HOFFSET(event_t, length), H5T_NATIVE_FLOAT);
    H5Tinsert(memtype, "mean", HOFFSET(event_t, mean), H5T_NATIVE_FLOAT);
    H5Tinsert(memtype, "stdv", HOFFSET(event_t, stdv), H5T_NATIVE_FLOAT);
    H5Tinsert(memtype, "pos", HOFFSET(event_t, pos), H5T_NATIVE_INT);

    // File representation
    hid_t filetype = H5Tcreate(H5T_COMPOUND, 4 * 5);
    if (filetype < 0) {
        warnx("Failed to create file representation for event table %s:%d.",
              __FILE__, __LINE__);
        goto clean2;
    }

    H5Tinsert(filetype, "start", 0, H5T_STD_U64LE);
    H5Tinsert(filetype, "length", 4, H5T_IEEE_F32LE);
    H5Tinsert(filetype, "mean", 4 * 2, H5T_IEEE_F32LE);
    H5Tinsert(filetype, "stdv", 4 * 3, H5T_IEEE_F32LE);
    H5Tinsert(filetype, "pos", 4 * 4, H5T_STD_I32LE);

    // Create dataset
    const hsize_t dims = et.n;
    hid_t space = H5Screate_simple(1, &dims, NULL);
    if (space < 0) {
        warnx("Failed to allocate dataspace for event table %s:%d.", __FILE__,
              __LINE__);
        goto clean3;
    }
    // Enable compression if available
    hid_t properties = H5P_DEFAULT;
    if (compression_level > 0) {
        properties = H5Pcreate(H5P_DATASET_CREATE);
        if (properties < 0) {
            warnx
                ("Failed to create properties structure to write out compressed data structure.\n");
            properties = H5P_DEFAULT;
        } else {
            H5Pset_shuffle(properties);
            H5Pset_deflate(properties, compression_level);
            H5Pset_chunk(properties, 1, &chunk_size);
        }
    }

    hid_t dset =
        H5Dcreate(hdf5file, readname, filetype, space, H5P_DEFAULT, properties,
                  H5P_DEFAULT);
    if (dset < 0) {
        warnx("Failed to create dataset for event table %s:%d.", __FILE__,
              __LINE__);
        goto clean4;
    }
    // Write data
    herr_t writeret =
        H5Dwrite(dset, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, et.event);
    if (writeret < 0) {
        warnx("Failed to write dataset for event table %s:%d.", __FILE__,
              __LINE__);
    }

 clean4:
    if (H5P_DEFAULT != properties) {
        H5Pclose(properties);
    }
    H5Dclose(dset);
 clean3:
    H5Sclose(space);
 clean2:
    H5Tclose(filetype);
 clean1:
    H5Tclose(memtype);
}

void write_annotated_raw(hid_t hdf5file, const char *readname,
                         const raw_table rt, hsize_t chunk_size,
                         int compression_level) {
    return;
}

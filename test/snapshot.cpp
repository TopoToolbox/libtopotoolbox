#include <gdal.h>
#undef NDEBUG

#include <gdal_priv.h>

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace tt {
extern "C" {
#include "topotoolbox.h"
}
}  // namespace tt

#include "profiler.h"
Profiler prof;

template <typename T, GDALDataType S>
void load_data_from_file(const std::string& filename, std::vector<T>& output,
                         float& cellsize, std::array<ptrdiff_t, 2>& dims) {
  GDALDataset* dataset = GDALDataset::Open(filename.data(), GA_ReadOnly);

  GDALRasterBand* poBand = dataset->GetRasterBand(1);
  dims[0] = poBand->GetXSize();
  dims[1] = poBand->GetYSize();

  output.resize(dims[0] * dims[1]);

  assert(poBand->RasterIO(GF_Read, 0, 0, dims[0], dims[1], output.data(),
                          dims[0], dims[1], S, 0, 0) == CE_None);

  double adfTransform[6];
  dataset->GetGeoTransform(adfTransform);
  cellsize = adfTransform[1];

  GDALClose(dataset);
}

template <typename T, GDALDataType S>
void write_data_to_file(const std::string& dst_filename,
                        const std::string& src_filename, std::vector<T>& output,
                        std::array<ptrdiff_t, 2>& dims) {
  GDALDataset* src_ds = GDALDataset::Open(src_filename.data(), GA_ReadOnly);
  assert(src_ds);

  GDALDriver* driver = src_ds->GetDriver();
  GDALDataset* dst_ds =
      driver->Create(dst_filename.data(), dims[0], dims[1], 1, S, NULL);
  assert(dst_ds);

  const OGRSpatialReference* ref = src_ds->GetSpatialRef();
  assert(dst_ds->SetSpatialRef(ref) == CE_None);

  double gt[6];
  assert(src_ds->GetGeoTransform(gt) == CE_None);
  assert(dst_ds->SetGeoTransform(gt) == CE_None);

  GDALRasterBand* poBand = dst_ds->GetRasterBand(1);
  assert(poBand->RasterIO(GF_Write, 0, 0, dims[0], dims[1], output.data(),
                          dims[0], dims[1], S, 0, 0) == CE_None);

  GDALClose(dst_ds);
  GDALClose(src_ds);
}

struct SnapshotData {
  std::filesystem::path path;

  std::array<ptrdiff_t, 2> dims;
  float cellsize;

  std::vector<float> dem;
  std::vector<float> filled_dem;
  std::vector<int32_t> flats;
  std::vector<int32_t> sills;
  std::vector<float> hs;
  std::vector<float> dilate3x3;
  std::vector<float> dilate_diag;
  std::vector<float> erode3x3;
  std::vector<float> erode_diag;

  // Output arrays
  std::vector<float> test_dem;
  std::vector<float> test_filled_dem;
  std::vector<int32_t> test_flats;
  std::vector<int32_t> test_sills;
  std::vector<float> test_hs;
  std::vector<float> test_nx;
  std::vector<float> test_ny;
  std::vector<float> test_nz;
  std::vector<float> test_dilate3x3;
  std::vector<float> test_dilate_diag;
  std::vector<float> test_erode3x3;
  std::vector<float> test_erode_diag;

  // Intermediate arrays
  std::vector<uint8_t> bc;
  std::vector<float> value_filter_tmp;

  SnapshotData(const std::filesystem::path& snapshot_path) {
    path = snapshot_path;

    GDALAllRegister();

    if (exists(snapshot_path / "dem.tif")) {
      load_data_from_file<float, GDT_Float32>(snapshot_path / "dem.tif", dem,
                                              cellsize, dims);
      test_dem = dem;  // Create a copy of the DEM to modify in fillsinks
      assert(test_dem.size() == dem.size());
    }

    std::array<ptrdiff_t, 2> dims_check;
    float cellcheck;
    if (exists(snapshot_path / "dilate3x3.tif")) {
      load_data_from_file<float, GDT_Float32>(snapshot_path / "dilate3x3.tif",
                                              dilate3x3, cellcheck, dims_check);
      assert(dims[0] == dims_check[0] && dims[1] == dims_check[1]);
    }

    if (exists(snapshot_path / "dilate_diag.tif")) {
      load_data_from_file<float, GDT_Float32>(snapshot_path / "dilate_diag.tif",
                                              dilate_diag, cellcheck,
                                              dims_check);
      assert(dims[0] == dims_check[0] && dims[1] == dims_check[1]);
    }

    if (exists(snapshot_path / "erode3x3.tif")) {
      load_data_from_file<float, GDT_Float32>(snapshot_path / "erode3x3.tif",
                                              erode3x3, cellcheck, dims_check);
      assert(dims[0] == dims_check[0] && dims[1] == dims_check[1]);
    }

    if (exists(snapshot_path / "erode_diag.tif")) {
      load_data_from_file<float, GDT_Float32>(
          snapshot_path / "erode_diag.tif", erode_diag, cellcheck, dims_check);
      assert(dims[0] == dims_check[0] && dims[1] == dims_check[1]);
    }

    if (exists(snapshot_path / "fillsinks.tif")) {
      load_data_from_file<float, GDT_Float32>(
          snapshot_path / "fillsinks.tif", filled_dem, cellcheck, dims_check);
      assert(dims[0] == dims_check[0] && dims[1] == dims_check[1]);
    }

    if (exists(snapshot_path / "identifyflats_flats.tif")) {
      load_data_from_file<int32_t, GDT_Int32>(
          snapshot_path / "identifyflats_flats.tif", flats, cellcheck,
          dims_check);
      assert(dims[0] == dims_check[0] && dims[1] == dims_check[1]);
    }

    if (exists(snapshot_path / "identifyflats_sills.tif")) {
      load_data_from_file<int32_t, GDT_Int32>(
          snapshot_path / "identifyflats_sills.tif", sills, cellcheck,
          dims_check);
      assert(dims[0] == dims_check[0] && dims[1] == dims_check[1]);
    }

    if (exists(snapshot_path / "hillshade.tif")) {
      load_data_from_file<float, GDT_Float32>(snapshot_path / "hillshade.tif",
                                              hs, cellcheck, dims_check);
      assert(dims[0] == dims_check[0] && dims[1] == dims_check[1]);
    }

    // Allocate and resize output and intermediate arrays
    if (dem.size() > 0) {
      if (erode3x3.size() > 0) {
        test_erode3x3.resize(dims[0] * dims[1]);
        value_filter_tmp.resize(dims[0] * dims[1]);
      }

      if (erode_diag.size() > 0) {
        test_erode_diag.resize(dims[0] * dims[1]);
      }

      if (dilate3x3.size() > 0) {
        test_dilate3x3.resize(dims[0] * dims[1]);
        value_filter_tmp.resize(dims[0] * dims[1]);
      }

      if (dilate_diag.size() > 0) {
        test_dilate_diag.resize(dims[0] * dims[1]);
      }

      if (filled_dem.size() > 0) {
        // Boundary conditions
        bc.resize(dims[0] * dims[1]);
        test_filled_dem.resize(dims[0] * dims[1]);
      }
      if (flats.size() > 0 && sills.size() > 0) {
        test_flats.resize(dims[0] * dims[1]);
      }

      if (hs.size() > 0) {
        test_nx.resize(dims[0] * dims[1]);
        test_ny.resize(dims[0] * dims[1]);
        test_nz.resize(dims[0] * dims[1]);
        test_hs.resize(dims[0] * dims[1]);
      }
    }
  }

  int test_min_filter_3x3() {
    ptrdiff_t se_dims[3] = {3, 3, 1};
    uint8_t se_full[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};

    tt::min_filter(test_erode3x3.data(), dem.data(), se_full, dims.data(),
                   se_dims);

    for (ptrdiff_t row = 0; row < dims[1]; row++) {
      for (ptrdiff_t col = 0; col < dims[0]; col++) {
        ptrdiff_t index = col + row * dims[0];
        if ((std::isnan(test_erode3x3[index]) != std::isnan(erode3x3[index])) ||
            // second isnan check not needed but added for completeness
            (!std::isnan(test_erode3x3[index]) &&
             !std::isnan(erode3x3[index]) &&
             test_erode3x3[index] != erode3x3[index])) {
          write_data_to_file<float, GDT_Float32>(path / "test_erode3x3.tif",
                                                 path / "erode3x3.tif",
                                                 test_erode3x3, dims);
          return -1;
        }
      }
    }
    return 0;
  }

  int test_min_filter_square_3x3() {
    uint8_t se_width = 3;

    tt::min_filter_square(test_erode3x3.data(), dem.data(),
                          value_filter_tmp.data(), se_width, dims.data());

    for (ptrdiff_t row = 0; row < dims[1]; row++) {
      for (ptrdiff_t col = 0; col < dims[0]; col++) {
        ptrdiff_t index = col + row * dims[0];
        if ((std::isnan(test_erode3x3[index]) != std::isnan(erode3x3[index])) ||
            // second isnan check not needed but added for completeness
            (!std::isnan(test_erode3x3[index]) &&
             !std::isnan(erode3x3[index]) &&
             test_erode3x3[index] != erode3x3[index])) {
          write_data_to_file<float, GDT_Float32>(
              path / "test_erode3x3_square.tif", path / "erode3x3.tif",
              test_erode3x3, dims);
          return -1;
        }
      }
    }
    return 0;
  }

  int test_min_filter_diag() {
    ptrdiff_t se_dims[3] = {3, 3, 1};
    uint8_t se_diag[9] = {1, 0, 1, 0, 1, 0, 1, 0, 1};

    tt::min_filter(test_erode_diag.data(), dem.data(), se_diag, dims.data(),
                   se_dims);

    for (ptrdiff_t row = 0; row < dims[1]; row++) {
      for (ptrdiff_t col = 0; col < dims[0]; col++) {
        ptrdiff_t index = col + row * dims[0];
        if ((std::isnan(test_erode_diag[index]) !=
             std::isnan(erode_diag[index])) ||
            // second isnan check not needed but added for completeness
            (!std::isnan(test_erode_diag[index]) &&
             !std::isnan(erode_diag[index]) &&
             test_erode_diag[index] != erode_diag[index])) {
          write_data_to_file<float, GDT_Float32>(path / "test_erode_diag.tif",
                                                 path / "erode_diag.tif",
                                                 test_erode_diag, dims);
          return -1;
        }
      }
    }
    return 0;
  }

  int test_max_filter_3x3() {
    ptrdiff_t se_dims[3] = {3, 3, 1};
    uint8_t se_full[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};

    tt::max_filter(test_dilate3x3.data(), dem.data(), se_full, dims.data(),
                   se_dims);

    for (ptrdiff_t row = 0; row < dims[1]; row++) {
      for (ptrdiff_t col = 0; col < dims[0]; col++) {
        ptrdiff_t index = col + row * dims[0];
        if ((std::isnan(test_dilate3x3[index]) !=
             std::isnan(dilate3x3[index])) ||
            // second isnan check not needed but added for completeness
            (!std::isnan(test_dilate3x3[index]) &&
             !std::isnan(dilate3x3[index]) &&
             test_dilate3x3[index] != dilate3x3[index])) {
          write_data_to_file<float, GDT_Float32>(path / "test_dilate3x3.tif",
                                                 path / "dilate3x3.tif",
                                                 test_dilate3x3, dims);
          return -1;
        }
      }
    }
    return 0;
  }

  int test_max_filter_square_3x3() {
    uint8_t se_width = 3;

    tt::max_filter_square(test_dilate3x3.data(), dem.data(),
                          value_filter_tmp.data(), se_width, dims.data());

    for (ptrdiff_t row = 0; row < dims[1]; row++) {
      for (ptrdiff_t col = 0; col < dims[0]; col++) {
        ptrdiff_t index = col + row * dims[0];
        if ((std::isnan(test_dilate3x3[index]) !=
             std::isnan(dilate3x3[index])) ||
            // second isnan check not needed but added for completeness
            (!std::isnan(test_dilate3x3[index]) &&
             !std::isnan(dilate3x3[index]) &&
             test_dilate3x3[index] != dilate3x3[index])) {
          write_data_to_file<float, GDT_Float32>(
              path / "test_dilate3x3_square.tif", path / "dilate3x3.tif",
              test_dilate3x3, dims);
          return -1;
        }
      }
    }
    return 0;
  }

  int test_max_filter_diag() {
    ptrdiff_t se_dims[3] = {3, 3, 1};
    uint8_t se_diag[9] = {1, 0, 1, 0, 1, 0, 1, 0, 1};

    tt::max_filter(test_dilate_diag.data(), dem.data(), se_diag, dims.data(),
                   se_dims);

    for (ptrdiff_t row = 0; row < dims[1]; row++) {
      for (ptrdiff_t col = 0; col < dims[0]; col++) {
        ptrdiff_t index = col + row * dims[0];
        if ((std::isnan(test_dilate_diag[index]) !=
             std::isnan(dilate_diag[index])) ||
            // second isnan check not needed but added for completeness
            (!std::isnan(test_dilate_diag[index]) &&
             !std::isnan(dilate_diag[index]) &&
             test_dilate_diag[index] != dilate_diag[index])) {
          write_data_to_file<float, GDT_Float32>(path / "test_dilate_diag.tif",
                                                 path / "dilate_diag.tif",
                                                 test_dilate_diag, dims);
          return -1;
        }
      }
    }
    return 0;
  }

  int test_fillsinks() {
    // Initialize bcs
    for (ptrdiff_t j = 0; j < dims[1]; j++) {
      for (ptrdiff_t i = 0; i < dims[0]; i++) {
        if (isnan(dem[j * dims[0] + i])) {
          bc[j * dims[0] + i] = 1;
          test_dem[j * dims[0] + i] = -INFINITY;
        } else {
          if (i == 0 || i == dims[0] - 1 || j == 0 || j == dims[1] - 1) {
            // 1 on the boundaries
            bc[j * dims[0] + i] = 1;
          } else {
            // 0 on the interior
            bc[j * dims[0] + i] = 0;
          }
        }
      }
    }

    tt::fillsinks(test_filled_dem.data(), test_dem.data(), bc.data(),
                  dims.data());

    for (ptrdiff_t j = 0; j < dims[1]; j++) {
      for (ptrdiff_t i = 0; i < dims[0]; i++) {
        if (!isnan(filled_dem[j * dims[0] + i])) {
          if (test_filled_dem[j * dims[0] + i] != filled_dem[j * dims[0] + i]) {
            write_data_to_file<float, GDT_Float32>(path / "test_fillsinks.tif",
                                                   path / "fillsinks.tif",
                                                   test_filled_dem, dims);
            return -1;
          }
        }
      }
    }
    return 0;
  }

  int test_identifyflats() {
    // identifyflats
    //
    // Use the snapshot filled DEM rather than the generated one in
    // case fillsinks fails.
    tt::identifyflats(test_flats.data(), filled_dem.data(), dims.data());

    for (ptrdiff_t j = 0; j < dims[1]; j++) {
      for (ptrdiff_t i = 0; i < dims[0]; i++) {
        if (flats[j * dims[0] + i] == 1 && !(test_flats[j * dims[0] + i] & 1)) {
          write_data_to_file<int32_t, GDT_Int32>(
              path / "test_identifyflats.tif", path / "identifyflats_flats.tif",
              test_flats, dims);
          return -1;
        }

        if (sills[j * dims[0] + i] == 1 && !(test_flats[j * dims[0] + i] & 2)) {
          write_data_to_file<int32_t, GDT_Int32>(
              path / "test_identifyflats.tif", path / "identifyflats_sills.tif",
              test_flats, dims);
          return -1;
        }
      }
    }
    return 0;
  }

  int test_hillshade() {
    // Azimuth and altitude are 315 and 60 degrees in radians
    // tt::hillshade requires azimuth to be in radians from the first
    // dimension towards the second dimension (a right-handed
    // coordinate system). The DEM data is loaded with the east
    // coordinates increasing in the first dimensions, and the north
    // coordinates decreasing in the second dimension. A bearing of
    // 315 degrees corresponds to 225 degrees (3.927 radians) in the
    // grid coordinate system.
    tt::hillshade(test_hs.data(), test_nx.data(), test_ny.data(), dem.data(),
                  3.9269908169872414, 1.047197551196598, cellsize, dims.data());

    for (ptrdiff_t j = 0; j < dims[1]; j++) {
      for (ptrdiff_t i = 0; i < dims[0]; i++) {
        float H = test_hs[j * dims[0] + i];

        if (fabsf(H - hs[j * dims[0] + i]) > 1e-4 ||
            (isnan(hs[j * dims[0] + i]) && !isnan(H))) {
          write_data_to_file<float, GDT_Float32>(path / "test_hillshade.tif",
                                                 path / "hillshade.tif",
                                                 test_hs, dims);
          return -1;
        }

        // Zero the hillshade array so we can reuse it for the
        // low-memory test
        test_hs[j * dims[0] + i] = 0.0;
      }
    }

    // Test the low-memory version of hillshade
    tt::hillshade_fused(test_hs.data(), dem.data(), 3.9269908169872414,
                        1.047197551196598, cellsize, dims.data());

    for (ptrdiff_t j = 0; j < dims[1]; j++) {
      for (ptrdiff_t i = 0; i < dims[0]; i++) {
        float H = test_hs[j * dims[0] + i];

        if (fabsf(H - hs[j * dims[0] + i]) > 1e-4) {
          write_data_to_file<float, GDT_Float32>(
              path / "test_hillshade_lowmem.tif", path / "hillshade.tif",
              test_hs, dims);
          return -1;
        }
      }
    }

    return 0;
  }

  int runtests() {
    std::cout << "    1..9" << std::endl;

    int result = 0;
    if (erode3x3.size() > 0) {
      if (test_min_filter_3x3() < 0) {
        result = -1;
        std::cout << "    not ok 1 - erode3x3" << std::endl;
      } else {
        std::cout << "    ok 1 - erode3x3" << std::endl;
      }
    }

    if (erode3x3.size() > 0) {
      if (test_min_filter_square_3x3() < 0) {
        result = -1;
        std::cout << "    not ok 2 - erode3x3_square" << std::endl;
      } else {
        std::cout << "    ok 2 - erode3x3_square" << std::endl;
      }
    }

    if (erode_diag.size() > 0) {
      if (test_min_filter_diag() < 0) {
        result = -1;
        std::cout << "    not ok 3 - erode_diag" << std::endl;
      } else {
        std::cout << "    ok 3 - erode_diag" << std::endl;
      }
    }

    if (dilate3x3.size() > 0) {
      if (test_max_filter_3x3() < 0) {
        result = -1;
        std::cout << "    not ok 4 - dilate3x3" << std::endl;
      } else {
        std::cout << "    ok 4 - dilate3x3" << std::endl;
      }
    }

    if (dilate3x3.size() > 0) {
      if (test_max_filter_square_3x3() < 0) {
        result = -1;
        std::cout << "    not ok 5 - dilate3x3_square" << std::endl;
      } else {
        std::cout << "    ok 5 - dilate3x3_square" << std::endl;
      }
    }

    if (dilate_diag.size() > 0) {
      if (test_max_filter_diag() < 0) {
        result = -1;
        std::cout << "    not ok 6 - dilate_diag" << std::endl;
      } else {
        std::cout << "    ok 6 - dilate_diag" << std::endl;
      }
    }

    // fillsinks
    if (test_filled_dem.size() > 0) {
      if (test_fillsinks() < 0) {
        result = -1;

        std::cout << "    not ok 7 - fillsinks" << std::endl;
      } else {
        std::cout << "    ok 7 - fillsinks" << std::endl;
      }
    }

    if (flats.size() > 0 && sills.size() > 0) {
      if (test_identifyflats() < 0) {
        result = -1;

        std::cout << "    not ok 8 - identifyflats" << std::endl;
      } else {
        std::cout << "    ok 8 - identifyflats" << std::endl;
      }
    }

    if (hs.size() > 0) {
      if (test_hillshade() < 0) {
        result = -1;
        std::cout << "    not ok 9 - hillshade" << std::endl;
      } else {
        std::cout << "    ok 9 - hillshade" << std::endl;
      }
    }

    return result;
  }
};

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Usage: snapshot <snapshot_directory>" << std::endl;
    return 0;
  }

  // Get the snapshot directory path from the command line arguments
  std::filesystem::path snapshot_dirpath(argv[1]);

  try {
    std::filesystem::directory_iterator dir(snapshot_dirpath);
  } catch (const std::filesystem::filesystem_error& e) {
    // Do not fail the test if the snapshots/data directory is not
    // present
    std::cout << "Snapshot directory does not exist" << std::endl;
    return 0;
  }

  std::ptrdiff_t test_count =
      std::distance(std::filesystem::directory_iterator{snapshot_dirpath},
                    std::filesystem::directory_iterator{});

  std::cout << "TAP version 14" << std::endl;
  std::cout << "1.." << test_count << std::endl;

  int result = 0;
  int test_id = 0;
  for (const auto& entry :
       std::filesystem::directory_iterator(snapshot_dirpath)) {
    SnapshotData data(entry.path());

    test_id++;
    int res = data.runtests();
    if (res < 0) {
      std::cout << "not ";
    }
    std::cout << "ok " << test_id;
    std::cout << " - " << entry.path().filename() << std::endl;

    result |= res;
  }
  return result;
}

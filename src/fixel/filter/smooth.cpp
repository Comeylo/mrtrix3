/*
 * Copyright (c) 2008-2018 the MRtrix3 contributors.
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, you can obtain one at http://mozilla.org/MPL/2.0/
 *
 * MRtrix3 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * For more details, see http://www.mrtrix.org/
 */


#include "fixel/filter/smooth.h"

#include "image_helpers.h"
#include "thread_queue.h"
#include "transform.h"
#include "fixel/helpers.h"

namespace MR
{
  namespace Fixel
  {
    namespace Filter
    {



      // TODO Add mask
      Smooth::Smooth (Image<Fixel::index_type> index_image,
                      const std::string& matrix_path,
                      const float fwhm,
                      const float threshold)
      {
        Fixel::check_index_image (index_image);
        const index_type num_fixels = Fixel::get_number_of_fixels (index_image);

        // Constants related to derivation of the smoothing matrix
        const float stdev = fwhm / 2.3548;
        const float gaussian_const1 = 1.0 / (stdev *  std::sqrt (2.0 * Math::pi));
        const float gaussian_const2 = -1.0 / (2.0 * stdev * stdev);

        // For generating the smoothing matrix, we need to be able to quickly
        //   calculate the distance between any pair of fixels
        vector<Eigen::Vector3> fixel_positions;
        fixel_positions.resize (num_fixels);
        Transform image_transform (index_image);
        for (auto i = Loop (index_image, 0, 3) (index_image); i; ++i) {
          Eigen::Vector3 vox ((default_type)index_image.index(0),
                              (default_type)index_image.index(1),
                              (default_type)index_image.index(2));
          vox = image_transform.voxel2scanner * vox;
          index_image.index(3) = 0;
          const index_type count = index_image.value();
          index_image.index(3) = 1;
          const index_type offset = index_image.value();
          for (index_type fixel_index = 0; fixel_index != count; ++fixel_index)
            fixel_positions[offset + fixel_index] = vox;
        }

        // We load the full connectivity matrix one line at a time,
        //   adding the spatial kernel to the smoothing
        matrix.reset (new Fixel::Matrix::norm_matrix_type (num_fixels, Fixel::Matrix::NormFixel()));

        class Source
        {
          public:
            Source (const std::string& filepath) :
                stream (filepath),
                progress ("Generating fixel data smoothing matrix"),
                index (0) { }
            Source (const Source&) = delete;
            bool operator() (std::pair<Fixel::index_type, std::string>& out)
            {
              out.first = index;
              if (std::getline (stream, out.second)) {
                ++index;
                ++progress;
                return true;
              }
              out.second.clear();
              return false;
            }
          private:
            std::ifstream stream;
            ProgressBar progress;
            Fixel::index_type index;
        };

        auto functor = [&] (const std::pair<Fixel::index_type, std::string>& in)
        {
          assert ((*matrix)[in.first].empty());
          const Fixel::Matrix::NormFixel input_fixel (Fixel::Matrix::parse_line<Fixel::Matrix::NormFixel> (in.second));
          Fixel::Matrix::NormFixel output_fixel;
          const Eigen::Vector3& pos (fixel_positions[in.first]);

          Fixel::Matrix::connectivity_value_type sum_weights = Fixel::Matrix::connectivity_value_type(0);
          for (const auto& it : input_fixel) {
            const default_type sq_distance = (fixel_positions[it.index()] - pos).squaredNorm();
            const Fixel::Matrix::connectivity_value_type weight = it.value() * gaussian_const1 * std::exp (gaussian_const2 * sq_distance);
            if (weight >= threshold) {
              output_fixel.push_back (Fixel::Matrix::NormElement (it.index(), weight));
              sum_weights += weight;
            }
          }

          if (sum_weights) {
            // Normalise smoothing weights
            const Fixel::Matrix::connectivity_value_type norm_factor = Fixel::Matrix::connectivity_value_type(1) / sum_weights;
            for (auto i : output_fixel)
              i.normalise (norm_factor);
          } else {
            // For a matrix intended for smoothing, want to give a fixel that is within the mask
            //   full self-connectivity even if it's not visited by any streamlines
            output_fixel.push_back (Fixel::Matrix::NormElement (in.first, Fixel::Matrix::connectivity_value_type (1)));
          }

          (*matrix)[in.first] = std::move (output_fixel);
          return true;
        };

        Source source (matrix_path);
        Thread::run_queue (source,
                           Thread::batch (std::pair<Fixel::index_type, std::string>()),
                           Thread::multi (functor));

      }




      void Smooth::operator() (Image<float>& input, Image<float>& output) const
      {
        Fixel::check_data_file (input);
        Fixel::check_data_file (output);

        check_dimensions (input, output);

        if (size_t (input.size(0)) != matrix->size())
          throw Exception ("Size of fixel data file \"" + input.name() + "\" (" + str(input.size(0)) +
                           ") does not match fixel connectivity matrix (" + str(matrix->size()) + ")");

        // Technically need to loop over axis 1; could be more than 1 parameter in the file
        for (auto l = Loop(1) (input, output); l; ++l) {
          for (size_t fixel = 0; fixel != matrix->size(); ++fixel) {
            input.index(0) = output.index(0) = fixel;
            if (std::isfinite (input.value())) {
              default_type value = 0.0, sum_weights = 0.0;
              for (const auto& it : (*matrix)[fixel]) {
                input.index(0) = it.index();
                if (std::isfinite (input.value())) {
                  value += input.value() * it.value();
                  sum_weights += it.value();
                }
              }
              if (sum_weights)
                output.value() = value / sum_weights;
              else
                output.value() = NaN;
            } else {
              output.value() = NaN;
            }
          }
        }
        input.index(1) = output.index(1) = 0;
      }



    }
  }
}


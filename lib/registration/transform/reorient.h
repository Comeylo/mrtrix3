/*
    Copyright 2013 Brain Research Institute, Melbourne, Australia

    Written by David Raffelt 2013

    This file is part of MRtrix.

    MRtrix is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MRtrix is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef __registration_transform_reorient_h__
#define __registration_transform_reorient_h__

#include "algo/threaded_loop.h"
#include "math/SH.h"
#include "math/least_squares.h"

namespace MR
{
  namespace Registration
  {
    namespace Transform
    {

      FORCE_INLINE Eigen::MatrixXd aPSF_weights_to_FOD_transform (const int num_SH, const Eigen::MatrixXd& directions) {
        Eigen::MatrixXd aPSF_matrix (num_SH, directions.cols());
        Math::SH::aPSF<float> aPSF_generator (Math::SH::LforN (num_SH));
        Eigen::Matrix<default_type, Eigen::Dynamic,1> aPSF;
        for (ssize_t i = 0; i < directions.cols(); ++i)
          aPSF_matrix.col(i) = aPSF_generator (aPSF, directions.col(i).head(3));
        return aPSF_matrix;
      }


      FORCE_INLINE Eigen::MatrixXd compute_reorient_transform (const ssize_t n_SH,
                                                               const transform_type& transform,
                                                               const Eigen::MatrixXd& directions) {
        return aPSF_weights_to_FOD_transform (n_SH, transform.linear().inverse() * directions)
               * Math::pinv(aPSF_weights_to_FOD_transform (n_SH, directions));
      }


      template <class FODImageType>
      class LinearKernel {

      public:
        LinearKernel (const Eigen::MatrixXd& transform) : transform (transform.cast <typename FODImageType::value_type> ()) {}


        void operator() (FODImageType& image) {
          image.index(3) = 0;
          if (image.value() > 0)  // only reorient voxels that contain a FODs
            image.row(3) = transform * image.row(3); // Eigen automatically takes care of the temporary
        }
         const Eigen::Matrix<typename FODImageType::value_type, Eigen::Dynamic, Eigen::Dynamic> transform;
      };


      template <class FODImageType>
      void reorient (FODImageType& fod_image,
                     const transform_type& transform,
                     const Eigen::MatrixXd& directions)
      {
        assert (directions.cols() > directions.rows());
        ThreadedLoop (fod_image, 0, 3)
            .run (LinearKernel<FODImageType>(compute_reorient_transform (fod_image.size(3), transform, directions)), fod_image);
      }


      template <class FODImageType>
      void reorient (const std::string progress_message,
                     FODImageType& fod_image,
                     const transform_type& transform,
                     const Eigen::MatrixXd& directions)
      {
        assert (directions.cols() > directions.rows());
        ThreadedLoop (progress_message, fod_image, 0, 3)
            .run (LinearKernel<FODImageType>(compute_reorient_transform (fod_image.size(3), transform, directions)), fod_image);
      }

    }
  }
}

#endif
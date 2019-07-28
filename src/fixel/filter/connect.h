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


#ifndef __fixel_filter_connect_h__
#define __fixel_filter_connect_h__

#include "image.h"
#include "fixel/matrix.h"
#include "fixel/filter/base.h"

#define DEFAULT_FIXEL_CONNECT_VALUE_THRESHOLD 0.5
#define DEFAULT_FIXEL_CONNECT_CONNECTIVITY_THRESHOLD 0.1

namespace MR
{
  namespace Fixel
  {
    namespace Filter
    {



      /** \addtogroup Filters
      @{ */

      /*! Perform a connected-component analysis of a fixel mask
       *
       * Typical usage:
       * \code
       * auto input = Image<float>::open (argument[0]);
       * Fixel::Matrix::norm_matrix_type matrix = Fixel::Matrix::load<Fixel::Matrix::norm_matrix_type> (argument[1]);
       * Fixel::Filter::Connect connect_filter (matrix);
       * auto output = Image::create<float> (argument[2], input);
       * smooth_filter (input, output);
       *
       * \endcode
       */

      class Connect : public Base
      { MEMALIGN (Connect)

        public:
          Connect (std::shared_ptr<Fixel::Matrix::norm_matrix_type> matrix,
                   const float value_threshold = DEFAULT_FIXEL_CONNECT_VALUE_THRESHOLD,
                   const float connectivity_threshold = DEFAULT_FIXEL_CONNECT_CONNECTIVITY_THRESHOLD) :
              matrix (matrix),
              value_threshold (value_threshold),
              connectivity_threshold (connectivity_threshold) { }

          void operator() (Image<float>& input, Image<float>& output) const;
          void set_value_threshold (const float value) { value_threshold = value; }
          void set_connectivity_threshold (const float value) { connectivity_threshold = value; }

        protected:
          std::shared_ptr<Fixel::Matrix::norm_matrix_type> matrix;
          float value_threshold, connectivity_threshold;
      };
    //! @}



    }
  }
}


#endif

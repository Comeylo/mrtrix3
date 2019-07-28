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


#ifndef __fixel_filter_smooth_h__
#define __fixel_filter_smooth_h__

#include "fixel/matrix.h"
#include "fixel/filter/base.h"

#define DEFAULT_FIXEL_SMOOTHING_FWHM 10.0
#define DEFAULT_FIXEL_SMOOTHING_THRESHOLD 0.01

namespace MR
{
  namespace Fixel
  {
    namespace Filter
    {



      /** \addtogroup Filters
      @{ */

      /*! Smooth fixel data using a combination of fixel-fixel connectivity and spatial distance.
       *
       * Typical usage:
       * \code
       * auto input = Image<float>::open (argument[0]);
       * auto index_image = Fixel::find_index_header (input.name()).get_image<Fixel::index_type>();
       * Fixel::Filter::Smooth smooth_filter (index_image, argument[1]);
       * auto output = Image::create<float> (argument[2], input);
       * smooth_filter (input, output);
       *
       * \endcode
       */

      class Smooth : public Base
      { MEMALIGN (Smooth)

        public:
          Smooth (Image<Fixel::index_type> index_image,
                  const std::string& matrix_path,
                  const float fwhm = DEFAULT_FIXEL_SMOOTHING_FWHM,
                  const float threshold = DEFAULT_FIXEL_SMOOTHING_THRESHOLD);

          void operator() (Image<float>& input, Image<float>& output) const;

        protected:
          std::shared_ptr<Fixel::Matrix::norm_matrix_type> matrix;
      };
    //! @}



    }
  }
}


#endif

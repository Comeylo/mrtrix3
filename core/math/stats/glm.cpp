/* Copyright (c) 2008-2017 the MRtrix3 contributors
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, you can obtain one at http://mozilla.org/MPL/2.0/.
 *
 * MRtrix is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * For more details, see http://www.mrtrix.org/.
 */


#include "math/stats/glm.h"

#include "debug.h"
#include "misc/bitset.h"
#include "thread_queue.h"

#define GLM_BATCH_SIZE 1024

namespace MR
{
  namespace Math
  {
    namespace Stats
    {
      namespace GLM
      {



        const char* const column_ones_description =
            "In some software packages, a column of ones is automatically added to the "
            "GLM design matrix; the purpose of this column is to estimate the \"global "
            "intercept\", which is the predicted value of the observed variable if all "
            "explanatory variables were to be zero. However there are rare situations "
            "where including such a column would not be appropriate for a particular "
            "experimental design. Hence, in MRtrix3 statistical inference commands, "
            "it is up to the user to determine whether or not this column of ones should "
            "be included in their design matrix, and add it explicitly if necessary. "
            "The contrast matrix must also reflect the presence of this additional column.";




        matrix_type solve_betas (const matrix_type& measurements, const matrix_type& design)
        {
          return design.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(measurements.transpose());
        }



        vector_type abs_effect_size (const matrix_type& measurements, const matrix_type& design, const Contrast& contrast)
        {
          if (contrast.is_F())
            return vector_type::Constant (measurements.rows(), NaN);
          else
            return matrix_type(contrast) * solve_betas (measurements, design);
        }

        matrix_type abs_effect_size (const matrix_type& measurements, const matrix_type& design, const vector<Contrast>& contrasts)
        {
          matrix_type result (measurements.rows(), contrasts.size());
          for (size_t ic = 0; ic != contrasts.size(); ++ic)
            result.col (ic) = abs_effect_size (measurements, design, contrasts[ic]);
          return result;
        }



        vector_type stdev (const matrix_type& measurements, const matrix_type& design)
        {
          matrix_type residuals = measurements.transpose() - design * solve_betas (measurements, design);
          residuals = residuals.array().pow (2.0);
          matrix_type one_over_dof (1, measurements.cols());
          one_over_dof.fill (1.0 / value_type(design.rows()-Math::rank (design)));
          return (one_over_dof * residuals).array().sqrt();
        }



        vector_type std_effect_size (const matrix_type& measurements, const matrix_type& design, const Contrast& contrast)
        {
          if (contrast.is_F())
            return vector_type::Constant (measurements.rows(), NaN);
          else
            return abs_effect_size (measurements, design, contrast).array() / stdev (measurements, design).array();
        }

        matrix_type std_effect_size (const matrix_type& measurements, const matrix_type& design, const vector<Contrast>& contrasts)
        {
          const auto stdev_reciprocal = vector_type::Ones (measurements.rows()).array() / stdev (measurements, design).array();
          matrix_type result (measurements.rows(), contrasts.size());
          for (size_t ic = 0; ic != contrasts.size(); ++ic)
            result.col (ic) = abs_effect_size (measurements, design, contrasts[ic]) * stdev_reciprocal;
          return result;
        }



        void all_stats (const matrix_type& measurements,
                        const matrix_type& design,
                        const vector<Contrast>& contrasts,
                        matrix_type& betas,
                        matrix_type& abs_effect_size,
                        matrix_type& std_effect_size,
                        vector_type& stdev)
        {
          ProgressBar progress ("calculating basic properties of default permutation");
          betas = solve_betas (measurements, design); ++progress;
          std::cerr << "Betas: " << betas.rows() << " x " << betas.cols() << ", max " << betas.array().maxCoeff() << "\n";
          abs_effect_size.resize (measurements.rows(), contrasts.size());
          for (size_t ic = 0; ic != contrasts.size(); ++ic) {
            if (contrasts[ic].is_F()) {
              abs_effect_size.col (ic).fill (NaN);
            } else {
              abs_effect_size.col (ic) = (matrix_type (contrasts[ic]) * betas).row (0);
            }
          }
          ++progress;
          std::cerr << "abs_effect_size: " << abs_effect_size.rows() << " x " << abs_effect_size.cols() << ", max " << abs_effect_size.array().maxCoeff() << "\n";
          matrix_type residuals = measurements.transpose() - design * betas;
          residuals = residuals.array().pow (2.0); ++progress;
          std::cerr << "residuals: " << residuals.rows() << " x " << residuals.cols() << ", max " << residuals.array().maxCoeff() << "\n";
          matrix_type one_over_dof (1, measurements.cols());
          one_over_dof.fill (1.0 / value_type(design.rows()-Math::rank (design)));
          std::cerr << "one_over_dof: " << one_over_dof.rows() << " x " << one_over_dof.cols() << ", max " << one_over_dof.array().maxCoeff() << "\n";
          VAR (design.rows());
          VAR (Math::rank (design));
          stdev = (one_over_dof * residuals).array().sqrt().row(0); ++progress;
          std::cerr << "stdev: " << stdev.size() << ", max " << stdev.array().maxCoeff() << "\n";
          std_effect_size = abs_effect_size.array() / stdev.array();  ++progress;
          std::cerr << "std_effect_size: " << std_effect_size.rows() << " x " << std_effect_size.cols() << ", max " << std_effect_size.array().maxCoeff() << "\n";
        }



        void all_stats (const matrix_type& measurements,
                        const matrix_type& fixed_design,
                        const vector<CohortDataImport>& extra_columns,
                        const vector<Contrast>& contrasts,
                        matrix_type& betas,
                        matrix_type& abs_effect_size,
                        matrix_type& std_effect_size,
                        vector_type& stdev)
        {
          if (extra_columns.empty()) {
            all_stats (measurements, fixed_design, contrasts, betas, abs_effect_size, std_effect_size, stdev);
            return;
          }

          class Source
          { NOMEMALIGN
            public:
              Source (const size_t num_elements) :
                  num_elements (num_elements),
                  counter (0),
                  progress (new ProgressBar ("calculating basic properties of default permutation", num_elements)) { }
              bool operator() (size_t& element_index)
              {
                element_index = counter++;
                if (counter >= num_elements) {
                  progress.reset();
                  return false;
                }
                assert (progress);
                ++(*progress);
                return true;
              }
            private:
              const size_t num_elements;
              size_t counter;
              std::unique_ptr<ProgressBar> progress;
          };

          class Functor
          { MEMALIGN(Functor)
            public:
              Functor (const matrix_type& data, const matrix_type& design_fixed, const vector<CohortDataImport>& extra_columns, const vector<Contrast>& contrasts,
                       matrix_type& betas, matrix_type& abs_effect_size, matrix_type& std_effect_size, vector_type& stdev) :
                  data (data),
                  design_fixed (design_fixed),
                  extra_columns (extra_columns),
                  contrasts (contrasts),
                  global_betas (betas),
                  global_abs_effect_size (abs_effect_size),
                  global_std_effect_size (std_effect_size),
                  global_stdev (stdev)
              {
                assert (design_fixed.cols() + extra_columns.size() == contrasts[0].cols());
              }
              bool operator() (const size_t& element_index)
              {
                const matrix_type element_data = data.row (element_index);
                matrix_type element_design (design_fixed.rows(), design_fixed.cols() + extra_columns.size());
                element_design.leftCols (design_fixed.cols()) = design_fixed;
                // For each element-wise design matrix column,
                //   acquire the data for this particular element, without permutation
                for (size_t col = 0; col != extra_columns.size(); ++col)
                  element_design.col (design_fixed.cols() + col) = (extra_columns[col]) (element_index);
                Math::Stats::GLM::all_stats (element_data, element_design, contrasts,
                                             local_betas, local_abs_effect_size, local_std_effect_size, local_stdev);
                global_betas.col (element_index) = local_betas;
                global_abs_effect_size.col (element_index) = local_abs_effect_size.col(0);
                global_std_effect_size.col (element_index) = local_std_effect_size.col(0);
                global_stdev[element_index] = local_stdev[0];
                return true;
              }
            private:
              const matrix_type& data;
              const matrix_type& design_fixed;
              const vector<CohortDataImport>& extra_columns;
              const vector<Contrast>& contrasts;
              matrix_type& global_betas;
              matrix_type& global_abs_effect_size;
              matrix_type& global_std_effect_size;
              vector_type& global_stdev;
              matrix_type local_betas, local_abs_effect_size, local_std_effect_size;
              vector_type local_stdev;
          };

          Source source (measurements.rows());
          Functor functor (measurements, fixed_design, extra_columns, contrasts,
                           betas, abs_effect_size, std_effect_size, stdev);
          Thread::run_queue (source, Thread::batch (size_t()), Thread::multi (functor));
        }






        Contrast::Partition Contrast::operator() (const matrix_type& design) const
        {
          // For now, let's do the most basic partitioning possible:
          //   Split design matrix column-wise depending on whether entries in the contrast matrix are all zero
          // TODO Later, may include config variables / compiler flags to change model partitioning technique
          matrix_type X, Z;
          const size_t nonzero_column_count = c.colwise().any().count();
          X.resize (design.rows(), nonzero_column_count);
          Z.resize (design.rows(), design.cols() - nonzero_column_count);
          ssize_t ix = 0, iz = 0;
          for (ssize_t ic = 0; ic != c.cols(); ++ic) {
            if (c.col (ic).any())
              X.col (ix++) = design.col (ic);
            else
              Z.col (iz++) = design.col (ic);
          }
          return Partition (X, Z);
        }








        TestFixed::TestFixed (const matrix_type& measurements, const matrix_type& design, const vector<Contrast>& contrasts) :
            TestBase (measurements, design, contrasts),
            pinvM (Math::pinv (M)),
            Rm (matrix_type::Identity (num_subjects(), num_subjects()) - (M*pinvM))
        {
          assert (contrasts[0].cols() == design.cols());
          // When the design matrix is fixed, we can pre-calculate the model partitioning for each contrast
          for (const auto c : contrasts)
            partitions.emplace_back (c (design));
        }



        void TestFixed::operator() (const matrix_type& shuffling_matrix, matrix_type& output) const
        {
          assert (shuffling_matrix.rows() == num_subjects());
          if (!(size_t(output.rows()) == num_elements() && size_t(output.cols()) == num_outputs()))
            output.resize (num_elements(), num_outputs());

          matrix_type beta, betahat;
          vector_type F;

          // Implement Freedman-Lane for fixed design matrix case
          // Each contrast needs to be handled explicitly on its own
          for (size_t ic = 0; ic != c.size(); ++ic) {

            // First, we perform permutation of the input data
            // In Freedman-Lane, the initial 'effective' regression against the nuisance
            //   variables, and permutation of the data, are done in a single step
            //VAR (perm_matrix.rows());
            //VAR (perm_matrix.cols());
            //VAR (partitions[ic].Rz.rows());
            //VAR (partitions[ic].Rz.cols());
            //VAR (y.rows());
            //VAR (y.cols());
            auto PRz = shuffling_matrix * partitions[ic].Rz;
            //VAR (PRz.rows());
            //VAR (PRz.cols());
            // TODO Re-attempt performing this as a single matrix multiplication across all elements
            matrix_type Sy (y.rows(), y.cols());
            for (ssize_t ie = 0; ie != y.rows(); ++ie)
              Sy.row (ie) = PRz * y.row (ie).transpose();
            //VAR (Sy.rows());
            //VAR (Sy.cols());
            // TODO Change measurements matrix convention to store data for each subject in a row;
            //   means data across subjects for a particular element appear in a column, which is contiguous

            // Now, we regress this shuffled data against the full model
            //VAR (pinvM.rows());
            //VAR (pinvM.cols());
            beta.noalias() = pinvM * Sy.transpose();
            //VAR (beta.rows());
            //VAR (beta.cols());
            //VAR (matrix_type(c[ic]).rows());
            //VAR (matrix_type(c[ic]).cols());
            betahat = matrix_type(c[ic]) * beta;
            //VAR (betahat.rows());
            //VAR (betahat.cols());
            //VAR (partitions[ic].X.rows());
            //VAR (partitions[ic].X.cols());
            //VAR (Rm.rows());
            //VAR (Rm.cols());
            auto XtX = partitions[ic].X.transpose()*partitions[ic].X;
            //VAR (XtX.rows());
            //VAR (XtX.cols());
            const default_type one_over_dof = 1.0 / (num_subjects() - partitions[ic].rank_x - partitions[ic].rank_z);
            auto residuals = Rm*Sy.transpose();
            //VAR (residuals.rows());
            //VAR (residuals.cols());
            vector_type temp3 = residuals.colwise().squaredNorm();
            //VAR (temp3.rows());
            //VAR (temp3.cols());
            // FIXME This should be giving a vector, not a matrix
            //auto temp4 = betahat.transpose() * (XtX * betahat) / c[ic].rank();
            //VAR (temp4.rows());
            //VAR (temp4.cols());
            //std::cerr << temp4 << "\n";
            F.resize (y.rows());
            for (ssize_t ie = 0; ie != y.rows(); ++ie) {
              vector_type this_betahat = betahat.col (ie);
              //VAR (this_betahat.size());
              auto temp2 = this_betahat.matrix() * (XtX * this_betahat.matrix()) / c[ic].rank();
              assert (temp2.rows() == 1);
              assert (temp2.cols() == 1);
              F[ie] = temp2 (0, 0) / (one_over_dof * temp3[ie]);
            }
            // TODO Try to use broadcasting here; it doesn't like having colwise() as the RHS argument
            //F = (betahat.transpose().rowwise() * ((partitions[ic].X.transpose()*partitions[ic].X) * betahat.colwise()) / c[ic].rank()) /
            //    ((Rm*Sy.transpose()).colwise().squaredNorm() / (num_subjects() - partitions[ic].rank_x - partitions[ic].rank_z));
            //VAR (F.size());

            // Put the results into the output matrix, replacing NaNs with zeroes
            // TODO Check again to see if this new statistic produces NaNs when input data are all zeroes
            // We also need to convert F to t if necessary
            for (ssize_t iF = 0; iF != F.size(); ++iF) {
              if (!std::isfinite (F[iF])) {
                output (iF, ic) = value_type(0);
              } else if (c[ic].is_F()) {
                output (iF, ic) = F[iF];
              } else {
                assert (betahat.rows() == 1);
                output (iF, ic) = std::sqrt (F[iF]) * (betahat (0, iF) > 0 ? 1.0 : -1.0);
              }
            }

          }
        }











        TestVariable::TestVariable (const vector<CohortDataImport>& importers,
                                    const matrix_type& measurements,
                                    const matrix_type& design,
                                    const vector<Contrast>& contrasts,
                                    const bool nans_in_data,
                                    const bool nans_in_columns) :
            TestBase (measurements, design, contrasts),
            importers (importers),
            nans_in_data (nans_in_data),
            nans_in_columns (nans_in_columns)
        {
          // Make sure that the specified contrasts reflect the full design matrix (with additional
          //   data loaded)
          assert (contrasts[0].cols() == M.cols() + ssize_t(importers.size()));
        }



        void TestVariable::operator() (const matrix_type& shuffling_matrix, matrix_type& output) const
        {
          if (!(size_t(output.rows()) == num_elements() && size_t(output.cols()) == num_outputs()))
            output.resize (num_elements(), num_outputs());

          // Let's loop over elements first, then contrasts in the inner loop
          for (ssize_t element = 0; element != y.rows(); ++element) {

            // For each element (row in y), need to load the additional data for that element
            //   for all subjects in order to construct the design matrix
            // Would it be preferable to pre-calculate and store these per-element design matrices,
            //   rather than re-generating them each time? (More RAM, less CPU)
            // No, most of the time that subject data will be memory-mapped, so pre-loading (in
            //   addition to the duplication of the fixed design matrix contents) would hurt bad
            matrix_type extra_data (num_subjects(), importers.size());
            for (ssize_t col = 0; col != ssize_t(importers.size()); ++col)
              extra_data.col(col) = importers[col] (element);

            // What can we do here that's common across all contrasts?
            // - Import the element-wise data
            // - Identify rows to be excluded based on NaNs in the design matrix
            // - Identify rows to be excluded based on NaNs in the input data
            //
            // Note that this is going to have to operate slightly differently to
            //   how it used to be done, i.e. via the permutation labelling vector,
            //   if we are to support taking the shuffling matrix as input to this functor
            // I think the approach will have to be:
            //   - Both NaNs in design matrix and NaNs in input data need to be removed
            //     in order to perform the initial regression against nuisance variables
            //   - Can then remove the corresponding _columns_ of the permutation matrix?
            //     No, don't think it's removal of columns; think it's removal of any rows
            //     that contain non-zero values in those columns
            //
            BitSet element_mask (M.rows(), true);
            if (nans_in_data) {
              for (ssize_t row = 0; row != y.rows(); ++row) {
                if (!std::isfinite (y (row, element)))
                  element_mask[row] = false;
              }
            }
            if (nans_in_columns) {
              for (ssize_t row = 0; row != extra_data.rows(); ++row) {
                if (!extra_data.row (row).allFinite())
                  element_mask[row] = false;
              }
            }
            const size_t finite_count = element_mask.count();

            // Do we need to reduce the size of our matrices / vectors
            //   based on the presence of non-finite values?
            matrix_type Mfull_masked;
            matrix_type perm_matrix_masked;
            vector_type y_masked;
            if (finite_count == num_subjects()) {

              Mfull_masked.resize (num_subjects(), num_factors());
              Mfull_masked.block (0, 0, num_subjects(), M.cols()) = M;
              Mfull_masked.block (0, M.cols(), num_subjects(), extra_data.cols()) = extra_data;
              perm_matrix_masked = shuffling_matrix;
              y_masked = y.row (element);

            } else {

              Mfull_masked.resize (finite_count, num_factors());
              y_masked.resize (finite_count);
              BitSet perm_matrix_mask (num_subjects(), true);
              ssize_t out_index = 0;
              for (size_t in_index = 0; in_index != num_subjects(); ++in_index) {
                if (element_mask[in_index]) {
                  Mfull_masked.block (out_index, 0, 1, M.cols()) = M.row (in_index);
                  Mfull_masked.block (out_index, M.cols(), 1, extra_data.cols()) = extra_data.row (in_index);
                  y_masked[out_index] = y (element, in_index);
                  ++out_index;
                } else {
                  // Any row in the permutation matrix that contains a non-zero entry
                  //   in the column corresponding to in_row needs to be removed
                  //   from the permutation matrix
                  for (ssize_t perm_row = 0; perm_row != shuffling_matrix.rows(); ++perm_row) {
                    if (shuffling_matrix (perm_row, in_index))
                      perm_matrix_mask[perm_row] = false;
                  }
                }
              }
              assert (out_index == ssize_t(finite_count));
              assert (perm_matrix_mask.count() == ssize_t(finite_count));
              // Only after we've reduced the design matrix do we now reduce the permutation matrix
              perm_matrix_masked.resize (finite_count, num_subjects());
              out_index = 0;
              for (size_t in_index = 0; in_index != num_subjects(); ++in_index) {
                if (perm_matrix_mask[in_index])
                  perm_matrix_masked.row (out_index++) = shuffling_matrix.row (in_index);
              }
              assert (out_index == finite_count);

            }
            assert (Mfull_masked.allFinite());

            const matrix_type pinvMfull_masked = Math::pinv (Mfull_masked);

            const matrix_type Rm = matrix_type::Identity (finite_count, finite_count) - (Mfull_masked*pinvMfull_masked);

            matrix_type beta, betahat;
            vector_type F;

            // We now have our permutation (shuffling) matrix and design matrix prepared,
            //   and can commence regressing the partitioned model of each contrast
            for (size_t ic = 0; ic != c.size(); ++ic) {

              const auto partition = c[ic] (Mfull_masked);

              // Now that we have the individual contrast model partition for these data,
              //   the rest of this function should proceed similarly to the fixed
              //   design matrix case
              // TODO Consider functionalising the below; should be consistent between fixed and variable
              const matrix_type Sy = perm_matrix_masked * partition.Rz * y_masked.matrix();
              beta.noalias() = Sy * pinvMfull_masked;
              betahat.noalias() = beta * matrix_type(c[ic]);
              F = (betahat.transpose() * (partition.X.inverse()*partition.X) * betahat / c[ic].rank()) /
                  ((Rm*Sy).squaredNorm() / (finite_count - partition.rank_x - partition.rank_z));

              for (ssize_t iF = 0; iF != F.size(); ++iF) {
                if (!std::isfinite (F[iF])) {
                  output (iF, ic) = value_type(0);
                } else if (c[ic].is_F()) {
                  output (iF, ic) = F[iF];
                } else {
                  assert (betahat.cols() == 1);
                  output (iF, ic) = std::sqrt (F[iF]) * (betahat (iF, 0) > 0 ? 1.0 : -1.0);
                }
              }

            } // End looping over contrasts

          } // End looping over elements
        }



        matrix_type TestVariable::default_design (const size_t index) const
        {
          matrix_type output (M.rows(), M.cols() + importers.size());
          output.block (0, 0, M.rows(), M.cols()) = M;
          for (size_t i = 0; i != importers.size(); ++i)
            output.col (M.cols() + i) = importers[i] (index);
          return output;
        }



      }
    }
  }
}


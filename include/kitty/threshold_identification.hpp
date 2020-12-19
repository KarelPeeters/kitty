/* kitty: C++ truth table library
 * Copyright (C) 2017-2020  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file threshold_identification.hpp
  \brief Threshold logic function identification

  \author CS-472 2020 Fall students
*/

#pragma once

#include <vector>
#include <lpsolve/lp_lib.h>
#include "traits.hpp"
#include "isop.hpp"

namespace kitty
{

namespace util
{

//TODO start using this again
/// Returns 1 for pos unate, -1 for neg unate and 0 for binate
template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
int is_unate( const TT& tt, uint8_t var_index )
{
  auto const tt0 = cofactor0( tt, var_index );
  auto const tt1 = cofactor1( tt, var_index );

  if ( is_const0( tt0 & ~tt1 ) )
    return 1;
  if ( is_const0( tt1 & ~tt0 ) )
    return -1;
  return 0;
}

} // namespace util

/*! \brief Threshold logic function identification

  Given a truth table, this function determines whether it is a threshold logic function (TF)
  and finds a linear form if it is. A Boolean function is a TF if it can be expressed as

  f(x_1, ..., x_n) = \sum_{i=1}^n w_i x_i >= T

  where w_i are the weight values and T is the threshold value.
  The linear form of a TF is the vector [w_1, ..., w_n; T].

  \param tt The truth table
  \param plf Pointer to a vector that will hold a linear form of `tt` if it is a TF.
             The linear form has `tt.num_vars()` weight values and the threshold value
             in the end.
  \return `true` if `tt` is a TF; `false` if `tt` is a non-TF.
*/
template<typename TT, typename = std::enable_if_t<is_complete_truth_table<TT>::value>>
bool is_threshold( const TT& tt_orig, std::vector<int64_t>* plf = nullptr )
{
  auto num_vars = tt_orig.num_vars();

  //convert to positive unate
  auto tt_star( tt_orig );
  std::vector<bool> flipped;

  for ( uint8_t var_index = 0u; var_index < num_vars; var_index++ )
  {
    auto x = util::is_unate( tt_orig, var_index );

    if ( x == 0 )
    {
      return false;
    }
    else if ( x == -1 )
    {
      flip_inplace( tt_star, var_index );
    }
    flipped.push_back(x == -1);
  }

  //construct ILP
  lprec* lp = make_lp( 0, num_vars + 1 );
  set_col_name( lp, 1, "T" );
  for ( int i = 0; i < num_vars + 1; i++ )
  {
    set_int( lp, i + 1, true );
  }

  auto col_indices = new int[num_vars + 1];
  auto row_values = new double[num_vars + 1];

  //-T as first col
  col_indices[0] = 1;
  row_values[0] = -1;

  set_add_rowmode( lp, true );

  std::vector<cube> on_sets = isop( tt_star );
  std::vector<cube> off_sets = isop( ~tt_star );

  for ( auto on_cube : on_sets )
  {
    std::cout << "on cube:  ";
    on_cube.print( num_vars );
    std::cout << std::endl;

    int j = 1;

    for ( int i = 0; i < num_vars; i++ )
    {
      if ( on_cube.get_mask( i ) && on_cube.get_bit( i ) )
      {
        col_indices[j] = i + 2;
        row_values[j] = 1;
        j++;
      }
    }

    add_constraintex( lp, j, row_values, col_indices, GE, 0 );
  }

  for ( auto off_cube : off_sets )
  {
    int j = 1;

    std::cout << "off cube: ";
    off_cube.print( num_vars );
    std::cout << std::endl;

    for ( int i = 0; i < num_vars; i++ )
    {
      if ( !off_cube.get_mask( i ) && !off_cube.get_bit( i ) )
      {
        col_indices[j] = i + 2;
        row_values[j] = 1;
        j++;
      }
    }

    add_constraintex( lp, j, row_values, col_indices, LE, -1 );
  }

  set_add_rowmode( lp, false );

  //add objective
  set_minim( lp );
  for ( int i = 0; i < num_vars + 1; i++ )
  {
    col_indices[i] = i + 1;
    row_values[i] = 1;
  }
  set_obj_fnex( lp, num_vars + 1, row_values, col_indices );

  write_LP( lp, stdout );

  //solve LP
  auto result = solve( lp );
  bool solved = false;

  if ( result == OPTIMAL )
  {

    print_solution( lp, num_vars + 1 );

    solved = true;
    if ( plf )
    {
      //convert solution to expected return format
      plf->clear();
      get_variables( lp, row_values );
      auto t_result = row_values[0];
      for ( int i = 0; i < num_vars; i++ )
      {
        plf->push_back( (int64_t)row_values[i + 1] );
        if ( flipped[i] )
        {
          t_result -= plf->back();
          plf->back() *= -1;
        }
      }
      plf->push_back( (int64_t)t_result );
    }
  }

  delete[] col_indices;
  delete[] row_values;

  return solved;
}

} /* namespace kitty */

// $Id: MitalasIntegration.h,v 1.3.4.3 2008-08-21 00:56:19 manav Exp $

// C/ C++ includes
#include <cmath>

// FESystem includes
#include "FESystem/FESystemNumbers.h"
#include "FESystem/FESystemExceptions.h"

// libmesh includes
#include "geom/node.h"

#ifdef ENABLE_ACCELERATE_FRAMEWORK
// Apple accelerate framework
#include "vecLib.h"
#endif // ENABLE_ACCELERATE_FRAMEWORK

inline double MitalasEdgeFactor(const Point& line1_p1,
                                const Point& line1_q1,
                                const Point& line2_p2,
                                const Point& line2_q2)
{
  // calculate the lengths for this pair
  static Point unit_vec1, unit_vec2, u_point, du_point, temp_vec1,
  p1, p2, q1, q2;
  double length1, length2, integral, d1, d2, f1, f2, T, J, V, cos_lambda,
    cos_theta, omega, du, factor = 0.0;
  unsigned int sum  = 0;
  static const unsigned int n_divs = 10;
  static double cos_vec_in[2*n_divs], cos_vec_out[2*n_divs],
    log_vec_in[2*n_divs], log_vec_out[2*n_divs], V_vec[n_divs];
  
  bool orientation_switch = false;
  
  unit_vec1.zero(); unit_vec2.zero(); u_point.zero(); du_point.zero(); temp_vec1.zero();
  length1=0.0; length2=0.0; integral=0.0; 
  d1=0.0; d2=0.0; f1=0.0; f2=0.0;
  T=0.0; J=0.0; V=0.0; cos_lambda=0.0; cos_theta=0.0; omega=0.0; du=0.0;
  
  // calculate the lengths
  unit_vec1.assign(line1_q1); unit_vec1 -= line1_p1; length1 = unit_vec1.size();
  unit_vec2.assign(line2_q2); unit_vec2 -= line2_p2; length2 = unit_vec2.size();
  
  // calculate the unit vectors
  factor = 1.0/length1; unit_vec1 *= factor;
  factor = 1.0/length2; unit_vec2 *= factor;
  
  // if the two vectors have opposite orientation, switch the point of 
  // one of the vectors and then proceede with the calculations. Finally, multiply the 
  // integral with -1.
  p1.assign(line1_p1);
  q1.assign(line1_q1);
  if (unit_vec1 * unit_vec2 <= 0.0)
    {
    p2.assign(line2_q2);
    q2.assign(line2_p2);
    unit_vec2 *= -1.0;
    orientation_switch = true;
    }
  else
    {
    p2.assign(line2_p2);
    q2.assign(line2_q2);
    }
  
  // calculate the unit vecs to find if the points p1,q1 lie on the other line segment
  temp_vec1.assign(p1); temp_vec1 -= p2; d1 = temp_vec1.size();
  if (d1 > FESystemNumbers::Epsilon)
    {
    factor = 1.0/d1; temp_vec1 *= factor;
    if (fabs(fabs(temp_vec1 * unit_vec2) -1) <= FESystemNumbers::Epsilon)
      sum += 1;
    }
  else
    sum += 1; // since the point coincides

  temp_vec1.assign(q1); temp_vec1 -= p2; f1 = temp_vec1.size();
  if (f1 > FESystemNumbers::Epsilon)
    {
    factor = 1.0/f1; temp_vec1 *= factor;
    if (fabs(fabs(temp_vec1 * unit_vec2) -1) <= FESystemNumbers::Epsilon)
      sum += 2;
    }
  else
    sum += 2; // since the point coincides
    

  temp_vec1.assign(p1); temp_vec1 -= q2; d2 = temp_vec1.size();
  temp_vec1.assign(q1); temp_vec1 -= q2; f2 = temp_vec1.size();
  
  integral = 0.0;
  du = length2 / (1.0 * n_divs);
  // this is the starting point for intgration
  u_point.assign(p2);
  
  // this is the increment in the integration point
  // this is scaled to half the value in order to increment the initial
  // point. Later is it increased to the actual value through multiplication by 2.0
  du_point.assign(unit_vec2);
  du_point *= (du/2.0);
  
  // now increment the point to the middle of the first interval
  u_point += du_point;
  du_point *= 2.0;

  switch (sum)
    {
    case 0:
      {
        // calculate all terms by numerical integration
        for (unsigned int i=0; i < n_divs; i++)
          {
          // calculate T, J, V, cos lambda, cos theta, omega
          temp_vec1.assign(p1); temp_vec1 -= u_point;
          J = temp_vec1.size();
          factor = 1.0/J; temp_vec1 *= factor;
          cos_theta = (-unit_vec1 * temp_vec1);
          if (fabs(cos_theta) > 1.0)
            {
            if (cos_theta < 0.0)
              cos_theta = -1.0;
            else
              cos_theta = 1.0;
            }
          
          // calculate the point for calculation of V
          temp_vec1.assign(p1);
          temp_vec1 += (unit_vec1 * cos_theta * J);
          temp_vec1 -= u_point;
          V = temp_vec1.size();

          // repeat as above for T and cos_lambda
          temp_vec1.assign(q1); temp_vec1 -= u_point;
          T = temp_vec1.size();
          factor = 1.0/T; temp_vec1 *= factor;
          cos_lambda = (unit_vec1 * temp_vec1); 
          if (fabs(cos_lambda) > 1.0)
            {
            if (cos_lambda < 0.0)
              cos_lambda = -1.0;
            else
              cos_lambda = 1.0;
            }

          // cos theta starts at 0, and cos lambda starts at n_divs
          cos_vec_in[i] = cos_lambda;
          cos_vec_in[i+n_divs] = cos_theta;
          
          // T starts at 0, and J starts at n_divs
          log_vec_in[i] = T;
          log_vec_in[i+n_divs] = J;
          V_vec[i] = V;
          // increment the point 
          u_point += du_point;
          }
        

#ifdef ENABLE_ACCELERATE_FRAMEWORK
        // calculate the cos values
        int_val = 2*n_divs;
        vvacos(cos_vec_out, cos_vec_in, &int_val);
        
        // calculate the log values
        vvlog(log_vec_out, log_vec_in, &int_val);
#else 
        for (unsigned int i=0; i<2*n_divs; i++)
          {
          cos_vec_out[i] = acos(cos_vec_in[i]);
          log_vec_out[i] = log(log_vec_in[i]);
          }
#endif  // ENABLE_ACCELERATE_FRAMEWORK 

        // calculate the integral values
        for (unsigned int i=0; i<n_divs; i++)
          {
          //
          // the original statement which has been vectorized is
          // omega = FESystemNumbers::Pi - (acos(cos_lambda) + acos(cos_theta));
          //
          omega = FESystemNumbers::Pi - (cos_vec_out[i] + cos_vec_out[i+n_divs]);
          
          //
          // the original statement which has been vectorized is
          // integral += (T * cos_lambda * log(T) + J * cos_theta * log(J) + V * omega) * du;
          //
          integral += (log_vec_in[i] * cos_vec_in[i] * log_vec_out[i] +
                       log_vec_in[i+n_divs] * cos_vec_in[i+n_divs] * log_vec_out[i+n_divs] +
                       V_vec[i] * omega) * du;
          }
      }
      break;
      
    case 1:
      {
        // direct calculation of cos theta term
        cos_theta = unit_vec1 * unit_vec2;
        if (fabs(cos_theta) > 1.0)
          {
          if (cos_theta < 0.0)
            cos_theta = -1.0;
          else
            cos_theta = 1.0;
          }
        if (d2 > 0.0)
          integral += cos_theta * (d2*d2 * (0.5 * log(d2) - 0.25));
        if (d1 > 0.0)
          integral -= cos_theta * (d1*d1 * (0.5 * log(d1) - 0.25));
        
        // if the point p1 is on q2 or beyond it on unit_vec2, then
        // cos theta needs to be scaled by -1.0
        if ((d2 == 0.0 || d1 > d2) && 
            ((d1 + d2 - length2) >= 0.0))
          cos_theta *= -1.0;
        
        // if p1 does not lie between the two points p2 and q2, the perform the 
        // integration, else break the line segment
        if ((d1 + d2 - length2) >= 0.0)
          {
          // calculate all other terms by numerical integration
          for (unsigned int i=0; i < n_divs; i++)
            {
            // calculate T, V, cos lambda, omega
            temp_vec1.assign(q1); temp_vec1 -= u_point;
            T = temp_vec1.size();
            factor = 1.0/T; temp_vec1 *= factor;
            cos_lambda = (unit_vec1 * temp_vec1); 
            if (fabs(cos_lambda) > 1.0)
              {
              if (cos_lambda < 0.0)
                cos_lambda = -1.0;
              else
                cos_lambda = 1.0;
              }
            
            // calculate the point for calculation of V
            temp_vec1.assign(q1);
            temp_vec1 -= (unit_vec1 * cos_lambda * T);
            temp_vec1 -= u_point;
            V = temp_vec1.size();
            
            // cos theta starts at 0, and cos lambda starts at n_divs
            cos_vec_in[i] = cos_lambda;

            // T starts at 0, and J starts at n_divs
            log_vec_in[i] = T;
            V_vec[i] = V;

            // increment the point 
            u_point += du_point;
            }

          // cos_theta is constant, hence, needs to be evaluated only once
          cos_vec_in[n_divs] = cos_theta;

#ifdef ENABLE_ACCELERATE_FRAMEWORK
          // calculate the cos values
          int_val = 1+n_divs;
          vvacos(cos_vec_out, cos_vec_in, &int_val);
          
          // calculate the log values
          int_val = n_divs;
          vvlog(log_vec_out, log_vec_in, &int_val);
#else
          for (unsigned int i=0; i<n_divs; i++)
            {
            cos_vec_out[i] = acos(cos_vec_in[i]);
            log_vec_out[i] = log(log_vec_in[i]);
            }
          cos_vec_out[n_divs] = acos(cos_vec_in[n_divs]);
#endif // ENABLE_ACCELERATE_FRAMEWORK
          
          // calculate the integral values
          for (unsigned int i=0; i<n_divs; i++)
            {
            //
            // the original statement which has been vectorized is
            // omega = FESystemNumbers::Pi - (acos(cos_lambda) + acos(cos_theta));
            //
            omega = FESystemNumbers::Pi - (cos_vec_out[i] + cos_vec_out[n_divs]);
            
            //
            // the original statement which has been vectorized is
            // integral += (T * cos_lambda * log(T) + V * omega) * du;
            //
            integral += (log_vec_in[i] * cos_vec_in[i] * log_vec_out[i] +
                         V_vec[i] * omega) * du;
            }
          }
        else
          {
          // needs to be implemented
          AssertThrow(false, ExcInternalError());
          }
      }
      break;
      
    case 2:
      {
        // direct calculation of cos lambda term
        cos_lambda = unit_vec1 * unit_vec2;
        if (fabs(cos_lambda) > 1.0)
          {
          if (cos_lambda < 0.0)
            cos_lambda = -1.0;
          else
            cos_lambda = 1.0;
          }
        
        if (f2 > 0.0)
          integral -= cos_lambda * (f2*f2 * (0.5 * log(f2) - 0.25));
        if (f1 > 0.0)
          integral += cos_lambda * (f1*f1 * (0.5 * log(f1) - 0.25));
        
        // if the point p1 is on q2 or beyond it on unit_vec2, then
        // cos theta needs to be scaled by -1.0
        if ((f1 == 0.0 || f2 > f1) && 
            ((f1 + f2 - length2) >= 0.0))
          cos_lambda *= -1.0;
        
        // if p1 does not lie between the two points p2 and q2, the perform the 
        // integration, else break the line segment
        if ((f1 + f2 - length2) >= 0.0)
          {
          // calculate all terms by numerical integration
          for (unsigned int i=0; i < n_divs; i++)
            {
            // calculate J, V, cos theta, omega
            temp_vec1.assign(p1); temp_vec1 -= u_point;
            J = temp_vec1.size();
            factor = 1.0/J; temp_vec1 *= factor;
            cos_theta = (-unit_vec1 * temp_vec1);
            if (fabs(cos_theta) > 1.0)
              {
              if (cos_theta < 0.0)
                cos_theta = -1.0;
              else
                cos_theta = 1.0;
              }
            
            // calculate the point for calculation of V
            temp_vec1.assign(p1);
            temp_vec1 += (unit_vec1 * cos_theta * J);
            temp_vec1 -= u_point;
            V = temp_vec1.size();
            
            // cos theta starts at 0, and cos lambda starts at n_divs
            cos_vec_in[1+i] = cos_theta;
            
            // T starts at 0, and J starts at n_divs
            log_vec_in[i] = J;
            V_vec[i] = V;
            // increment the point 
            u_point += du_point;
            }

          // cos_lambda is constant, hence is evaluated only once
          cos_vec_in[0] = cos_lambda;
          
#ifdef ENABLE_ACCELERATE_FRAMEWORK
          // calculate the cos values
          int_val = 1+n_divs;
          vvacos(cos_vec_out, cos_vec_in, &int_val);
          
          // calculate the log values
          int_val = n_divs;
          vvlog(log_vec_out, log_vec_in, &int_val);
#else
          for (unsigned int i=0; i<n_divs; i++)
            {
            cos_vec_out[i] = acos(cos_vec_in[i]);
            log_vec_out[i] = log(log_vec_in[i]);
            }
          cos_vec_out[n_divs] = acos(cos_vec_in[n_divs]);
#endif // ENABLE_ACCELERATE_FRAMEWORK
          
          // calculate the integral values
          for (unsigned int i=0; i<n_divs; i++)
            {
            //
            // the original statement which has been vectorized is
            // omega = FESystemNumbers::Pi - (acos(cos_lambda) + acos(cos_theta));
            //
            omega = FESystemNumbers::Pi - (cos_vec_out[0] + cos_vec_out[1+i]);
            
            //
            // the original statement which has been vectorized is
            // integral += (J * cos_theta * log(J) + V * omega) * du;
            //
            integral += (log_vec_in[i] * cos_vec_in[1+i] * log_vec_out[i] +
                         V_vec[i] * omega) * du;
            }
          }
        else
          {
          AssertThrow(false, ExcInternalError());
          }
      }
      break;

    case 3:
      {
        // direct calculation of all terms
        if (d2 > 0.0)
          integral += (d2*d2 * (0.5 * log(d2) - 0.25));
        if (d1 > 0.0)
          integral -= (d1*d1 * (0.5 * log(d1) - 0.25));
        if (f2 > 0.0)
          integral -= (f2*f2 * (0.5 * log(f2) - 0.25));
        if (f1 > 0.0)
          integral += (f1*f1 * (0.5 * log(f1) - 0.25));
        
      }
      break;
      
    default:
      Assert(false, ExcInternalError());
    }
  
  integral -= length1 * length2;
  integral *= (unit_vec1 * unit_vec2);
  if (orientation_switch)
    integral *= -1.0;
  
#ifdef DEBUG
  if (integral != integral)
    {
    std::cout << line1_p1 << line1_q1 << std::endl;
    std::cout << line2_p2 << line2_q2 << std::endl;
    
    std::cout << "line integral : " << integral << std::endl << std::endl;
//    std::cout << MitalasEdgeFactor(line1_p1, line1_q1, line2_p2, line2_q2);
    }
#endif
  
  return integral;
}


inline double NumIntEdgeFactor(const Point& line1_p1,
                               const Point& line1_q1,
                               const Point& line2_p2,
                               const Point& line2_q2)
{
  // calculate the lengths for this pair
  static Point unit_vec1, unit_vec2, u_point, du_point, e_point, de_point_b2, 
  de_point, temp_vec1,
  p1, p2, q1, q2;
  double length1, length2, integral, d1, d2, f1, f2, T, J, V, cos_lambda,
    cos_theta, omega, du;
  static const unsigned int n_divs = 50;
  bool orientation_switch = false;
  
  unit_vec1.zero(); unit_vec2.zero(); u_point.zero(); du_point.zero(); 
  e_point.zero(); temp_vec1.zero(); de_point.zero();
  length1=0.0; length2=0.0; integral=0.0; 
  d1=0.0; d2=0.0; f1=0.0; f2=0.0;
  T=0.0; J=0.0; V=0.0; cos_lambda=0.0; cos_theta=0.0; omega=0.0; du=0.0;
  
  // calculate the lengths
  unit_vec1.assign(line1_q1); unit_vec1 -= line1_p1; length1 = unit_vec1.size();
  unit_vec2.assign(line2_q2); unit_vec2 -= line2_p2; length2 = unit_vec2.size();
  
  // calculate the unit vectors
  unit_vec1 /= length1;
  unit_vec2 /= length2;
  
  // if the two vectors have opposite orientation, switch the point of 
  // one of the vectors and then proceede with the calculations. Finally, multiply the 
  // integral with -1.
  p1.assign(line1_p1);
  q1.assign(line1_q1);
  if (unit_vec1 * unit_vec2 <= 0.0)
    {
    p2.assign(line2_q2);
    q2.assign(line2_p2);
    unit_vec2 *= -1.0;
    orientation_switch = true;
    }
  else
    {
    p2.assign(line2_p2);
    q2.assign(line2_q2);
    }
  
  integral = 0.0;
  du = length2 / (1.0 * n_divs);
  // this is the starting point for intgration
  u_point.assign(p2);
  
  // this is the increment in the integration point
  // this is scaled to half the value in order to increment the initial
  // point. Later is it increased to the actual value through multiplication by 2.0
  du_point.assign(unit_vec2);
  du_point *= (du/2.0);
  
  de_point.assign(unit_vec1);
  de_point *= (du);
  de_point_b2.assign(unit_vec1);
  de_point_b2 *= (du/2);
  
  // now increment the point to the middle of the first interval
  u_point += du_point;
  du_point *= 2.0;

  // calculate all terms by numerical integration
  for (unsigned int i=0; i < n_divs; i++)
    {

    e_point.assign(p1);
    e_point += de_point_b2;
    for (unsigned int i=0; i < n_divs; i++)
      {
      temp_vec1.assign(e_point);
      temp_vec1 -= u_point;
      
      integral += log(temp_vec1.size()) * du * du;
      
      e_point += de_point;
      }    
    // once all the integration is done, increment the point 
    u_point += du_point;
    }
  
  integral *= (unit_vec1 * unit_vec2);
  if (orientation_switch)
    integral *= -1.0;
  
//  std::cout << line1_p1 << line1_q1 << std::endl;
//  std::cout << line2_p2 << line2_q2 << std::endl;
//
//  std::cout << "line integral : " << integral << std::endl;
  
  return integral;
}



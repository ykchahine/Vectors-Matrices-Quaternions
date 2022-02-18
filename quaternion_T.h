//
//  quaternion_T.h
//  linear_algebra
//
//  Created by William McCarthy on 172//21.
//

#ifndef quaternion_T_h
#define quaternion_T_h

#include <cmath>
#include "vector3d_T.h"
#include "matrix3d_T.h"

template <typename T> class quaternion;
template <typename T> using quat = class quaternion<T>;
typedef quat<double> quatD;

template <typename T>
class quaternion {
public:
//w is the scalar
  quaternion(T w_=T(), T x_=T(), T y_=T(), T z_=T())
  : w(w_), x(x_), y(y_), z(z_) { }

  static quaternion i() { return (0, 1, 0, 0); }
  static quaternion j() { return (0, 0, 1, 0); }
  static quaternion k() { return (0, 0, 0, 1); }

  static double ii() {return -1;}
  static double jj() {return -1;}
  static double kk() {return -1;}
  static double ijk() {return -1;}

  static quaternion ij() {return k();}
  static quaternion jk() {return i();}
  static quaternion ki() {return j();}

  static quaternion ji() {return -k();}
  static quaternion kj() {return -i();}
  static quaternion ik() {return -j();}

  friend quaternion operator+(const quaternion& a, const quaternion& b){
     //w + i*x + j*y + k*z + 
     
     return quaternion( (a.w + b.w), (a.x + b.x), (a.y + b.y), (a.z + b.z) );

  }
  friend quaternion operator-(const quaternion& a, const quaternion& b){
     
     return quaternion( (a.w - b.w), (a.x - b.x), (a.y - b.y), (a.z - b.z) );

  }
  friend quaternion operator*(const quaternion& a, const quaternion& b) {
     return quaternion((a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z),
                              (a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y),
                              (a.w * b.y + a.y * b.w + a.z * b.x - a.x * b.z),
                              (a.w * b.z + a.z * b.w + a.x * b.y - a.y * b.w));
  }
  friend quaternion operator+(const quaternion& q, T k){
    return quaternion( (q.w + k), (q.x + k), (q.y + k), (q.z + k) );
  }
  friend quaternion operator+(T k, const quaternion& q) { return(q + k); }

  friend quaternion operator-(const quaternion& q, T k){
    return quaternion( (q.w - k), (q.x - k), (q.y - k), (q.z - k) );

  }
  friend quaternion operator-(T k, const quaternion& q) { return -(q-k); }

  friend quaternion operator*(const quaternion& q, T k){
    return quaternion( (q.w * k), (q.x * k), (q.y * k), (q.z * k) );
  } //4 * (1,2,3,4)
  friend quaternion operator*(T k, const quaternion& q) {return (q * k);}
  friend quaternion operator/(const quaternion& q, T k){
    while(k != 0) {
    return quaternion( (q.w / k), (q.x / k), (q.y / k), (q.z / k) );
  }
  
  }
  
   quaternion operator-() const {
    return quaternion(-w, -x, -y, -z);
  } //make everything negative

  friend bool operator==(const quaternion& q, const quaternion& r){
    return (q.w == r.w && q.x == r.x && q.y == r.y && q.z == r.z);
  }
  friend bool operator!=(const quaternion& q, const quaternion& r) {return ! (q == r);}
  
  vector3d<T> vector() const {return vector3d<T>("quat", 3, {x, y, z});}
  
  T scalar() {return w;}

  quaternion unit_scalar() const { return quaternion(1.0, vector()); }

  quaternion conjugate() const { return quaternion(w, -x, -y, -z); }

  quaternion inverse() { return conjugate() / (magnitude() * magnitude()); }

  quaternion unit() const { return quaternion() / magnitude(); }

  double norm() const { 
    quaternion<T> a = *this;
    return sqrt((a.w * a.w)  + (a.x * a.x) + (a.y * a.y) + (a.z * a.z));
  }
  double magnitude() const {
    quaternion<T> a = *this;
    return a.norm(); 
  }

  double dot(const quaternion& v) const{

    return w * v.w + vector().dot(v.vector());

  }

  double angle(const quaternion& v) const {

        quaternion<T> z = conjugate() * v;
        T zvnorm = z.vector().norm();
        T zscalar = z.scalar();
        double angle = atan2(zvnorm, zscalar);
        return angle * 180.0 / 3.1415;

  }

  matrix3d<T> rot_matrix() const{
    
     T w = w;
     T x = x;
     T y = y;
     T z = z;
     return matrix3d<T>("mat", 3, { -2*(y*y + z*z) + 1,  2*(x*y  - w*z), 2*(x*z  +  w*y),
                             2*(x*y  + w*z),      -2*(x*x + z*z) + 1,  2*(y*z  -  w*x),
                             2*(x*z  - w*y),       2*(y*z  + w*x),      -2*(x*x + y*y)+1) };
  }

 // rotates point pt (pt.x, pt.y, pt.z) about (axis.x, axis.y, axis.z) by theta
 static vec3 rotate(const vector3D& pt, const vector3D& axis, double theta);

 friend std::ostream& operator<<(std::ostream& os, const quaternion& q) {
   os << "Quat(";
   if (q ==  quaternion::i())  { return os <<  "i)"; }
   if (q == -quaternion::i())  { return os << "-i)"; }
   if (q ==  quaternion::j())  { return os <<  "j)"; }
   if (q == -quaternion::j())  { return os << "-j)"; }
   if (q ==  quaternion::k())  { return os <<  "k)"; }
   if (q == -quaternion::k())  { return os << "-k)"; }

   if (q.magnitude() == 0.0 && q.w == 0)   { return os << "0)"; }
   if (q.magnitude() == 0.0 && q.w == 0)   { return os << "0)"; }
   if (q.magnitude() == 1.0 && q.w == 1)   { return os << "1)"; }
   if (q.vector().magnitude() == 0.0)      { return os << q.w << ")"; }
   else { return os << q.w << q.vector() << ")"; }
 }

 static void run_tests();

private:
 T w, x, y, z;
};

void plane_rotation(const std::string& msg, const quatD& plane, const std::initializer_list<double>& li) {
 matrix3dD rotate = matrix3dD("rot_matrix", 3, li);
 assert(plane.rot_matrix() == rotate);
 std::cout << msg << " is: " << plane << plane.rot_matrix() << "\n";
}


std::string yes_or_no(bool condition) { return condition ? "YES" : "no"; }

template <typename T>
void quaternion<T>::run_tests() {
 quatD a(1, 2, 3, 4), b(4, 0, 0, 7), c(0, 1, 1, 0), d(0, 0, 1, 0);
 quatD e(0, 0, 0, 1), f(0, 0, 0, 0), g(1, 0, 0, 0), h(3, 0, 0, 0);
 std::cout << "a = " << a << ")\nb = " << b << ")\nc = " << c << ")\nd = " << d
           << ")\ne = " << e << ")\nf = " << f << ")\ng = " << g << ")\nh = " <<  h << "\n";

 std::cout << "c + d = " <<  c + d << "\nc + d + e = " << c + d + e << "\n";
 std::cout << "5 * h = " << 5 * h << "\nh * 5 = " << h * 5 << "\nh / 3.0 = " << h / 3.0 << "\n\n";

 std::cout << "h.magnitude() is " << h.magnitude() << "\nh.unit() is " << h.unit();
 std::cout << "g.unit() is " << g.unit() << "\na.unit() is " << a.unit() << ")\n\n";

 std::cout << "a.vector() is " << a.vector() << "\na.scalar() is " << a.scalar() << "\n";
 std::cout << "a.conjugate() is " << a.conjugate() << "\na.inverse() is " << a.inverse()
           << "\na * a.inverse() is " << a * a.inverse() << "\n\n";

 std::cout << "c == d is " << yes_or_no(c == d) << "\nc != d is " << yes_or_no(c != d);
 std::cout << "\ne == e is " << yes_or_no(e == e) << "\ne != e is " << yes_or_no(e != e) << "\n";

 std::cout << "\n\nquat.ij is: " << quatD::ij() << "\nquat.jk is: " << quatD::jk()
           << "\nquat.ki is: " << quatD::ki() << "\n";
 assert(quatD::ij() == quatD::k());
 assert(quatD::jk() == quatD::i());
 assert(quatD::ki() == quatD::j());

 std::cout << "\nquat.ji is: " << quatD::ji() << "\nquat.kj is: " << quatD::kj()
           << "\nquat.ik is: " << quatD::ik() << "\nquat.ijk is: " << quatD::ijk() << "\n";
 assert(quatD::ji() == -quatD::k());
 assert(quatD::kj() == -quatD::i());
 assert(quatD::ik() == -quatD::j());

 std::cout << "\nquat.ii is: " << quatD::ii() << "\nquat.jj is: " << quatD::jj()
           << "\nquat.kk is: " << quatD::kk() << "\n";
 assert(quatD::ii()  == -1);
 assert(quatD::jj()  == -1);
 assert(quatD::kk()  == -1);
 assert(quatD::ijk() == -1);

 std::cout << "\nangle (deg) between c and d is: " << c.angle(d) << "\n";
 quatD c_minus_d = c - d;
 std::cout << "c_minus_d is: " << c_minus_d;
 matrix3dD rot_matrix = c_minus_d.rot_matrix();
 std::cout << "rot_matrix of c_minus_d is: " << c_minus_d.rot_matrix() << "\n";

 double rad2_2 = sqrt(2)/2.0;
 std::cout << "// -------------- LEVEL FLIGHT -------------------')\n";
 plane_rotation("levelFlight(E)", quatD(1),                     {  1,  0,  0,   0,  1,  0,   0,  0,  1 });
 plane_rotation("levelFlight(N)", quatD(rad2_2, 0, rad2_2,  0), {  0,  0,  1,   0,  1,  0,  -1,  0,  0 });
 plane_rotation("levelFlight(W)", quatD(0,      0,  1,      0), { -1,  0,  0,   0,  1,  0,   0,  0, -1 });
 plane_rotation("levelFlight(S)", quatD(rad2_2, 0, -rad2_2, 0), {  0,  0, -1,   0,  1,  0,   1,  0,  0} );
 std::cout << "LEVEL FLIGHT assertions passed ..............................................\n";
 std::cout << "// --------- end LEVEL FLIGHT ------------------------)\n";

 std::cout << "// -------------- STRAIGHT UP -------------------')\n";
 plane_rotation("straightUp(E)", quatD(rad2_2, 0, 0, rad2_2),   {  0, -1,  0,   1,  0,  0,   0,  0,  1 } );
 plane_rotation("straightUp(N)", quatD(0.5, 0.5, 0.5, 0.5),     {  0,  0,  1,   1,  0,  0,   0,  1,  0 } );
 plane_rotation("straightUp(W)", quatD(0, rad2_2, rad2_2, 0),   {  0,  1,  0,   1,  0,  0,   0,  0, -1 } );
 plane_rotation("straightUp(S)", quatD(0.5, -0.5, -0.5, 0.5),   {  0,  0, -1,   1,  0,  0,   0, -1,  0 } );
 std::cout << "STRAIGHT UP assertions passed..............................................\n";
 std::cout << "// --------- end STRAIGHT UP ------------------------)\n\n";


 std::cout << "// -------------- STRAIGHT DOWN ------------------')\n";
 plane_rotation("straightDown(E)", quatD(rad2_2, 0, 0, -rad2_2), {  0,  1,  0,  -1,  0,  0,   0,  0,  1 } );
 plane_rotation("straightDown(E)", quatD(0.5, -0.5, 0.5, -0.5),  {  0,  0,  1,  -1,  0,  0,   0, -1,  0 });
 plane_rotation("straightDown(E)", quatD(0, -rad2_2, rad2_2, 0), {  0, -1,  0,  -1,  0,  0,   0,  0, -1 } );
 plane_rotation("straightDown(E)", quatD(0.5, 0.5, -0.5, -0.5),  {  0,  0, -1,  -1,  0,  0,   0,  1,  0 });
 std::cout << "STRAIGHT DOWN assertions passed..............................................\n";
 std::cout << "// --------- end STRAIGHT DOWN ----------------------)\n\n";


 std::cout << "\n\n -------- BANK/ROLL ----------------\n";
 std::cout << "\nBanking/Rolling 90 degrees left...\n";
 plane_rotation("plane_E_bankLeft90", quatD(rad2_2, rad2_2, 0, 0),  {  1,  0,  0,   0,  0, -1,   0,  1,  0 } );
 plane_rotation("plane_N_bankLeft90", quatD(0.5, 0.5, 0.5, -0.5),   {  0,  1,  0,   0,  0, -1,  -1,  0,  0 } );
 plane_rotation("plane_W_bankLeft90", quatD(0, 0, rad2_2, -rad2_2), { -1,  0,  0,   0,  0, -1,   0, -1,  0 }  );
 plane_rotation("plane_W_bankLeft90", quatD(0.5, 0.5, -0.5, 0.5),   {  0, -1,  0,   0,  0, -1,   1,  0,  0 } );
 std::cout << "ROLL 90 deg left assertions passed..............................................\n";

 std::cout << "\n\nBanking/Rolling 180 degrees...\n";
 plane_rotation("plane_E_bankLeft180", quatD(0, 1, 0, 0),            {  1,  0,  0,   0, -1,  0,   0,  0, -1 });
 plane_rotation("plane_N_bankLeft180", quatD(0, rad2_2, 0, -rad2_2), {  0,  0, -1,   0, -1,  0,  -1,  0,  0 });
 plane_rotation("plane_W_bankLeft180", quatD(0, 0, 0, 1),            { -1,  0,  0,   0, -1,  0,   0,  0,  1 });
 plane_rotation("plane_S_bankLeft180", quatD(0, rad2_2, 0, rad2_2),  {  0,  0,  1,   0, -1,  0,   1,  0,  0 });
 std::cout << "ROLL 180 degrees assertions passed..............................................\n";

 std::cout << "\n\nBanking/Rolling 90 degrees right...\n";
 plane_rotation("plane_E_bankRight90", quatD(rad2_2, -rad2_2, 0, 0), {  1,  0,  0,   0,  0,  1,   0, -1,  0 });
 plane_rotation("plane_N_bankRight90", quatD(0.5, -0.5, 0.5, 0.5),   {  0, -1,  0,   0,  0,  1,  -1,  0,  0 });
 plane_rotation("plane_W_bankRight90", quatD(0, 0, rad2_2, rad2_2),  { -1,  0,  0,   0,  0,  1,   0,  1,  0 });
 plane_rotation("plane_S_bankRight90", quatD(0.5, -0.5, -0.5, -0.5), {  0,  1,  0,   0,  0,  1,   1,  0,  0 });
 std::cout << "ROLL 90 deg right assertions passed..............................................\n";
 std::cout << "\n -------- end BANK/ROLL ----------------\n";

 std::cout << "\nALL PLANE ROTATION ASSERTIONS PASSED ............................................\n\n";

 std::cout << "SEE THIS WEBSITE for DETAILED DIAGRAMS on the TESTS of the PLANE's rotations\n";
 std::cout << "https://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToMatrix/examples/index.htm\n";
}


#endif /* quaternion_T_h */


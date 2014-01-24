#include <type_traits>

struct A{ typedef int I; I in;};
struct B{ typedef int I; I in;};

typedef A GridA;
typedef B GridB;

template <class T>
struct Info {

  bool ok(T& ) {
    return true;
  }
};

typedef Info<GridA> AI;
typedef Info<GridB> BI;


struct ProblemBase {

  virtual AI* info() = 0;
  virtual BI* subinfo() = 0;

};

struct Problem : public ProblemBase {

  virtual AI* info() { return new AI(); }
  virtual BI* subinfo() { return new BI(); }

};

Problem problem() { return Problem(); }

template <class J>
typename std::enable_if<std::is_same<GridA, J>::value, bool>::type is_neu(J& j) {
  return problem().info()->ok(j);
}

template <class J>
typename std::enable_if<!std::is_same<GridA, J>::value, bool>::type is_neu(J& j) {
  return problem().subinfo()->ok(j);
}

int main() {
  GridA a;
  GridB b;
  auto foo = is_neu(a);
  auto bar = is_neu(b);
}

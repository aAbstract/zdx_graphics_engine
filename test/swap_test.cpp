#include <iostream>

int main() {
  int a = 5;
  int b = 10;
  std::cout << "a -> " << a << ", b -> " << b << std::endl;
  std::swap(a, b);
  std::cout << "a -> " << a << ", b -> " << b << std::endl;
  return 0;
}
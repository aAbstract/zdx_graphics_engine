#include <iostream>

int get_item(int x, int y, int* arr, int rows, int cols) {
  if (x < rows && y < cols) {
    return arr[x * cols + y];
  } else {
    return -1;
  }
}

int main() {
  int rows = 3;
  int cols = 6;
  int* arr = new int[rows * cols];
  for (int i = 0; i < rows * cols; i++) {
    arr[i] = i * 10;
    std::cout << i << " -> " << i * 10 << std::endl;
  }
  std::cout << get_item(0, 2, arr, rows, cols) << std::endl;  // 20
  std::cout << get_item(0, 5, arr, rows, cols) << std::endl;  // 50
  std::cout << get_item(1, 3, arr, rows, cols) << std::endl;  // 90
  std::cout << get_item(1, 4, arr, rows, cols) << std::endl;  // 100
  std::cout << get_item(2, 3, arr, rows, cols) << std::endl;  // 150
  return 0;
}
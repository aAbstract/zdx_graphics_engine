#include <png.h>
#include <stdlib.h>
#include <iostream>

#define info_log std::cout << "[INFO]: "
#define err_log std::cout << "[ERROR]: "

struct pixel {
  uint8_t r, g, b;
};

class zdX_texture {
 private:
  pixel* raw_data;
  int32_t tex_width;
  int32_t tex_height;

 public:
  bool load_from_file(const char* file_name) {
    png_structp raw_png_data;
    png_infop raw_info_png_data;
    raw_png_data =
        png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!raw_png_data)
      return false;
    raw_info_png_data = png_create_info_struct(raw_png_data);
    if (!raw_info_png_data)
      return false;
    if (setjmp(png_jmpbuf(raw_png_data)))
      return false;
    FILE* f = fopen(file_name, "rb");
    if (f) {
      png_init_io(raw_png_data, f);
      png_read_info(raw_png_data, raw_info_png_data);
      png_byte color_type;
      png_byte bit_depth;
      png_bytep* row_pointers;
      this->tex_width = png_get_image_width(raw_png_data, raw_info_png_data);
      this->tex_height = png_get_image_height(raw_png_data, raw_info_png_data);
      color_type = png_get_color_type(raw_png_data, raw_info_png_data);
      bit_depth = png_get_bit_depth(raw_png_data, raw_info_png_data);
      if (bit_depth == 16)
        png_set_strip_16(raw_png_data);
      if (color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_palette_to_rgb(raw_png_data);
      if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
        png_set_expand_gray_1_2_4_to_8(raw_png_data);
      if (png_get_valid(raw_png_data, raw_info_png_data, PNG_INFO_tRNS))
        png_set_tRNS_to_alpha(raw_png_data);
      if (color_type == PNG_COLOR_TYPE_RGB ||
          color_type == PNG_COLOR_TYPE_GRAY ||
          color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_filler(raw_png_data, 0xFF, PNG_FILLER_AFTER);
      if (color_type == PNG_COLOR_TYPE_GRAY ||
          color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
        png_set_gray_to_rgb(raw_png_data);
      png_read_update_info(raw_png_data, raw_info_png_data);
      row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * this->tex_height);
      for (int y = 0; y < tex_height; y++) {
        row_pointers[y] = (png_byte*)malloc(
            png_get_rowbytes(raw_png_data, raw_info_png_data));
      }
      png_read_image(raw_png_data, row_pointers);
      this->raw_data = new pixel[this->tex_width * this->tex_height];
      for (int y = 0; y < this->tex_height; y++) {
        png_bytep row = row_pointers[y];
        for (int x = 0; x < this->tex_width; x++) {
          png_bytep px = &(row[x * 4]);
          this->raw_data[y * this->tex_width + x] = {px[0], px[1], px[2]};
        }
      }
      for (int y = 0; y < this->tex_height; y++)
        free(row_pointers[y]);
      free(row_pointers);
      png_destroy_read_struct(&raw_png_data, &raw_info_png_data, nullptr);
      fclose(f);
      return true;
    } else {
      return false;
    }
  }

  int get_width() const { return this->tex_width; }
  int get_height() const { return this->tex_height; }

  pixel get_pixel(int x, int y) const {
    if (x >= 0 && x < this->tex_width && y >= 0 && y < this->tex_height) {
      return this->raw_data[y * this->tex_width + x];
    } else {
      return {255, 0, 0};
    }
  }
};

int main() {
  return 0;
}
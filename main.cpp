#define OLC_PGE_APPLICATION

// #define DEBUG_MODE_RENDER
// #define DEBUG_MODE_PARAMS
// #define DEBUG_MODE_GFX
// #define WIREFRAME_MODE

#include "./olcPixelGameEngine.h"
#include "./zdx_engine_core.h"

#include <algorithm>
#include <vector>

#define info_log std::cout << "[INFO]: "
#define err_log std::cout << "[ERROR]: "

class zdx_engine : public olc::PixelGameEngine {
 private:
  // Fields
  float motion_speed = 25.0f;
  float rot_speed = 50.0f;
  float theta = 0.0f;
  int screen_width = 256;   // 256
  int screen_height = 200;  // 200

  zdx_core engine_core;
  _3d_surface _3d_model;
  zdx_texture ztex;

  void render(_3d_surface& surface) {
#ifdef DEBUG_MODE_RENDER
    info_log << "Rendering:-" << endl;
    for (int i = 0; i < surfce.tris.size(); i++) {
      std::cout << "T" << i << ": {";
      std::cout << surfce.tris[i].p1.x << ",";
      std::cout << surfce.tris[i].p1.y << ",";
      std::cout << surfce.tris[i].p1.z << ",";
      std::cout << "   ";
      std::cout << surfce.tris[i].p2.x << ",";
      std::cout << surfce.tris[i].p2.y << ",";
      std::cout << surfce.tris[i].p2.z << ",";
      std::cout << "   ";
      std::cout << surfce.tris[i].p3.x << ",";
      std::cout << surfce.tris[i].p3.y << ",";
      std::cout << surfce.tris[i].p3.z << ",";
      std::cout << " }\n";
    }
#endif
    if (surface.has_texture) {
      for (int i = 0; i < surface.tris.size(); i++) {
        std::vector<vec_info> temp;
        temp.reserve(3);
        temp.emplace_back(surface.tris[i].p1, surface.tris[i].tp1);
        temp.emplace_back(surface.tris[i].p2, surface.tris[i].tp2);
        temp.emplace_back(surface.tris[i].p3, surface.tris[i].tp3);
        this->draw_tex_tri(temp, this->ztex);

#ifdef WIREFRAME_MODE
        this->DrawTriangle(
            {(int)surface.tris[i].p1.x, (int)surface.tris[i].p1.y},
            {(int)surface.tris[i].p2.x, (int)surface.tris[i].p2.y},
            {(int)surface.tris[i].p3.x, (int)surface.tris[i].p3.y}, olc::WHITE);
#endif
      }
    } else {
      for (int i = 0; i < surface.tris.size(); i++) {
        this->FillTriangle(
            {(int)surface.tris[i].p1.x, (int)surface.tris[i].p1.y},
            {(int)surface.tris[i].p2.x, (int)surface.tris[i].p2.y},
            {(int)surface.tris[i].p3.x, (int)surface.tris[i].p3.y},
            surface.tris[i].face_color);

#ifdef WIREFRAME_MODE
        this->DrawTriangle(
            {(int)surface.tris[i].p1.x, (int)surface.tris[i].p1.y},
            {(int)surface.tris[i].p2.x, (int)surface.tris[i].p2.y},
            {(int)surface.tris[i].p3.x, (int)surface.tris[i].p3.y}, olc::WHITE);
#endif
      }
    }
  }

  void draw_tex_tri(std::vector<vec_info>& ps, const zdx_texture& tex_info) {
    std::sort(ps.begin(), ps.end(), [](const vec_info& a, const vec_info& b) {
      return a.p.y < b.p.y;
    });

    // First Half of The Triangle
    float delta_y1 = ps[1].p.y - ps[0].p.y;
    float delta_y2 = ps[2].p.y - ps[0].p.y;
    float delta_x1 = ps[1].p.x - ps[0].p.x;
    float delta_x2 = ps[2].p.x - ps[0].p.x;

    float delta_u1 = ps[1].tp.x - ps[0].tp.x;
    float delta_u2 = ps[2].tp.x - ps[0].tp.x;
    float delta_v1 = ps[1].tp.y - ps[0].tp.y;
    float delta_v2 = ps[2].tp.y - ps[0].tp.y;

    float dx_dy1 = 0, dx_dy2 = 0, du_dy1 = 0, du_dy2 = 0, dv_dy1 = 0,
          dv_dy2 = 0;

    if (delta_y2) {
      dx_dy2 = delta_x2 / delta_y2;
      du_dy2 = delta_u2 / delta_y2;
      dv_dy2 = delta_v2 / delta_y2;
    }
    if (delta_y1) {
      dx_dy1 = delta_x1 / delta_y1;
      du_dy1 = delta_u1 / delta_y1;
      dv_dy1 = delta_v1 / delta_y1;

      for (int i = 0; i <= (int)delta_y1; i++) {
        // x coordinate of the space start point
        int sxi = ps[0].p.x + i * dx_dy1;
        int sxf = ps[0].p.x + i * dx_dy2;

        // (x, y) coordinate of the texture start point
        float tui = ps[0].tp.x + i * du_dy1;
        float tvi = ps[0].tp.y + i * dv_dy1;

        // (x, y) coordinate of the texture end point
        float tuf = ps[0].tp.x + i * du_dy2;
        float tvf = ps[0].tp.y + i * dv_dy2;

        if (sxi > sxf) {
          std::swap(sxi, sxf);
          std::swap(tui, tuf);
          std::swap(tvi, tvf);
        }

        if ((sxf - sxi) != 0) {
          float du_dx = (tuf - tui) / (sxf - sxi);
          float dv_dx = (tvf - tvi) / (sxf - sxi);
          for (int j = 0; j <= (sxf - sxi); j++) {
            float tex_u = tui + j * du_dx;
            float tex_v = tvi + j * dv_dx;
            this->Draw(sxi + j, ps[0].p.y + i,
                       tex_info.get_pixel_norm(tex_u, tex_v));
          }
        }
      }
    }

    // Second Half of The Triangle
    float delta_y3 = ps[2].p.y - ps[1].p.y;
    float delta_x3 = ps[2].p.x - ps[1].p.x;

    float delta_u3 = ps[2].tp.x - ps[1].tp.x;
    float delta_v3 = ps[2].tp.y - ps[1].tp.y;

    if (delta_y3) {
      float dx_dy3 = delta_x3 / delta_y3;
      float du_dy3 = delta_u3 / delta_y3;
      float dv_dy3 = delta_v3 / delta_y3;

      for (int i = 0; i <= (int)delta_y3; i++) {
        // x coordinate of the space start point
        int sxi = ps[1].p.x + i * dx_dy3;
        int sxf = ps[0].p.x + (i + delta_y1) * dx_dy2;

        // (x, y) coordinate of the texture start point
        float tui = ps[1].tp.x + i * du_dy3;
        float tvi = ps[1].tp.y + i * dv_dy3;

        // (x, y) coordinate of the texture end point
        float tuf = ps[0].tp.x + (i + delta_y1) * du_dy2;
        float tvf = ps[0].tp.y + (i + delta_y1) * dv_dy2;

        if (sxi > sxf) {
          std::swap(sxi, sxf);
          std::swap(tui, tuf);
          std::swap(tvi, tvf);
        }

        if ((sxf - sxi) != 0) {
          float du_dx = (tuf - tui) / (sxf - sxi);
          float dv_dx = (tvf - tvi) / (sxf - sxi);
          for (int j = 0; j <= (sxf - sxi); j++) {
            float tex_u = tui + j * du_dx;
            float tex_v = tvi + j * dv_dx;
            this->Draw(sxi + j, ps[1].p.y + i,
                       tex_info.get_pixel_norm(tex_u, tex_v));
          }
        }
      }
    }
  }

 public:
  zdx_engine()
      : engine_core(
            60.0f,
            1000.0f,
            this->get_screen_width(),
            this->get_screen_height(),
            {{0.0f, 0.0f, 0.0f}, {0.0f, 0.0f, 1.0f}, {1.0f, 0.0f, 0.0f}},
            {{0.0f, 0.0f, 1.0f}},
            0.1f) {
    sAppName = "Camera Effect";
  }

  // Getters
  int get_screen_width() { return this->screen_width; }
  int get_screen_height() { return this->screen_height; }

  bool OnUserCreate() override {
    const char* model_file = "./spyro_tex.obj";
    const char* texture_file = "./Low.png";
    info_log << "Done Constructing Graphics Window\n";
    if (this->_3d_model.load_from_file(model_file, true) &&
        this->ztex.load_from_file(texture_file)) {
      return true;
    } else {
      err_log << "Asset Files Not Found\n";
      return false;
    }
  }

  bool OnUserUpdate(float fElapsedTime) override {
    // Update Frame Variables
    this->theta += fElapsedTime * this->rot_speed;
    if (this->theta >= 360) {
      this->theta = 0;
    }

    // Updata Engine Parameters
    if (this->GetKey(olc::UP).bHeld)
      this->engine_core.get_camera().cam_vec.y -=
          this->motion_speed * fElapsedTime;

    if (this->GetKey(olc::DOWN).bHeld)
      this->engine_core.get_camera().cam_vec.y +=
          this->motion_speed * fElapsedTime;

    if (this->GetKey(olc::RIGHT).bHeld)
      this->engine_core.yaw_cam(-1.0f * this->rot_speed * fElapsedTime);

    if (this->GetKey(olc::LEFT).bHeld)
      this->engine_core.yaw_cam(this->rot_speed * fElapsedTime);

    if (this->GetKey(olc::W).bHeld) {
      this->engine_core.get_camera().cam_vec.x +=
          this->motion_speed * fElapsedTime *
          this->engine_core.get_camera().z_vec.x;
      this->engine_core.get_camera().cam_vec.y +=
          this->motion_speed * fElapsedTime *
          this->engine_core.get_camera().z_vec.y;
      this->engine_core.get_camera().cam_vec.z +=
          this->motion_speed * fElapsedTime *
          this->engine_core.get_camera().z_vec.z;
    }

    if (this->GetKey(olc::S).bHeld) {
      this->engine_core.get_camera().cam_vec.x -=
          this->motion_speed * fElapsedTime *
          this->engine_core.get_camera().z_vec.x;
      this->engine_core.get_camera().cam_vec.y -=
          this->motion_speed * fElapsedTime *
          this->engine_core.get_camera().z_vec.y;
      this->engine_core.get_camera().cam_vec.z -=
          this->motion_speed * fElapsedTime *
          this->engine_core.get_camera().z_vec.z;
    }

    if (this->GetKey(olc::D).bHeld) {
      this->engine_core.get_camera().cam_vec.x +=
          this->motion_speed * fElapsedTime *
          this->engine_core.get_camera().x_vec.x;
      this->engine_core.get_camera().cam_vec.y +=
          this->motion_speed * fElapsedTime *
          this->engine_core.get_camera().x_vec.y;
      this->engine_core.get_camera().cam_vec.z +=
          this->motion_speed * fElapsedTime *
          this->engine_core.get_camera().x_vec.z;
    }

    if (this->GetKey(olc::A).bHeld) {
      this->engine_core.get_camera().cam_vec.x -=
          this->motion_speed * fElapsedTime *
          this->engine_core.get_camera().x_vec.x;
      this->engine_core.get_camera().cam_vec.y -=
          this->motion_speed * fElapsedTime *
          this->engine_core.get_camera().x_vec.y;
      this->engine_core.get_camera().cam_vec.z -=
          this->motion_speed * fElapsedTime *
          this->engine_core.get_camera().x_vec.z;
    }

    // Draw Frame
    // Clear Screen and Draw
    this->FillRect(0, 0, this->get_screen_width(), this->get_screen_height(),
                   olc::L);
    _3d_surface projected_model =
        this->engine_core.bulk_pipeline(this->_3d_model, this->theta);
    this->render(projected_model);

    // for (int i = 0; i < this->ztex.get_width(); i++) {
    //   for (int j = 0; j < this->ztex.get_height(); j++) {
    //     this->Draw(i, j, this->ztex.get_pixel(i, j));
    //   }
    // }
    return true;
  }
};

int main() {
  info_log << "ZDX 3D Engine\n";
  info_log << "Constructing Graphics Window...\n";
  zdx_engine engine;
  if (engine.Construct(engine.get_screen_width(), engine.get_screen_height(), 4,
                       4))
    engine.Start();

  return 0;
}
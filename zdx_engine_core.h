#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector>

// Engine Data Sturctures
struct vec_3d {
  float x, y, z, w = 1;
};

struct vec_2d {
  float x, y, w = 1;
};

struct mat_4x4 {
  float mat[4][4] = {0.0f};
};

struct triangle {
  vec_3d p1, p2, p3;
  vec_2d tp1, tp2, tp3;
  olc::Pixel face_color = 0;
};

struct vec_info {
  vec_3d p;
  vec_2d tp;
  vec_info(const vec_3d& _p, const vec_2d& _tp) : p(_p), tp(_tp) {}
};

struct camera {
  vec_3d cam_vec;
  vec_3d z_vec;
  vec_3d x_vec;
  vec_3d y_vec;
};

struct _3d_surface {
  std::vector<triangle> tris;
  bool has_texture;
  bool load_from_file(const char* file_name, bool _has_texture = false) {
    this->has_texture = _has_texture;
    std::ifstream file(file_name);
    if (!file.is_open())
      return false;
    std::vector<vec_3d> verts_pool;
    std::vector<vec_2d> tex_verts_pool;
    while (!file.eof()) {
      char line_buffer[64];
      file.getline(line_buffer, 64);
      std::stringstream reader;
      reader << line_buffer;
      if (!has_texture) {
        char d_type;
        reader >> d_type;
        if (d_type == 'v') {
          vec_3d temp_vec;
          reader >> temp_vec.x >> temp_vec.y >> temp_vec.z;
          verts_pool.push_back(temp_vec);
        } else if (d_type == 'f') {
          int vi[3];
          reader >> vi[0] >> vi[1] >> vi[2];
          triangle temp_tri;
          temp_tri.p1 = verts_pool[vi[0] - 1];
          temp_tri.p2 = verts_pool[vi[1] - 1];
          temp_tri.p3 = verts_pool[vi[2] - 1];
          this->tris.push_back(temp_tri);
        }
      } else {
        char d_type[2];
        reader >> d_type;
        if (d_type[0] == 'v' && d_type[1] == 't') {
          vec_2d temp_vec;
          reader >> temp_vec.x >> temp_vec.y;
          tex_verts_pool.push_back(temp_vec);
        } else if (d_type[0] == 'v') {
          vec_3d temp_vec;
          reader >> temp_vec.x >> temp_vec.y >> temp_vec.z;
          verts_pool.push_back(temp_vec);
        } else if (d_type[0] == 'f') {
          int vi[3];
          int tvi[3];
          std::string faces_params[3];
          reader >> faces_params[0] >> faces_params[1] >> faces_params[2];
          for (int i = 0; i < 3; i++) {
            std::stringstream temp1(faces_params[i]);
            std::string segment;
            std::vector<std::string> seglist;
            while (std::getline(temp1, segment, '/')) {
              seglist.push_back(segment);
            }
            if (seglist.size() == 0) {
              continue;
            }
            std::stringstream temp2(seglist[0]);
            int x;
            temp2 >> x;
            vi[i] = x;
            temp2 = std::stringstream(seglist[1]);
            temp2 >> x;
            tvi[i] = x;
          }
          triangle temp_tri;
          temp_tri.p1 = verts_pool[vi[0] - 1];
          temp_tri.p2 = verts_pool[vi[1] - 1];
          temp_tri.p3 = verts_pool[vi[2] - 1];
          temp_tri.tp1 = tex_verts_pool[tvi[0] - 1];
          temp_tri.tp2 = tex_verts_pool[tvi[1] - 1];
          temp_tri.tp3 = tex_verts_pool[tvi[2] - 1];
          this->tris.push_back(temp_tri);
        }
      }
    }
    return true;
  }
};

class zdx_core {
 private:
  // Transformation Matrices
  mat_4x4 proj_mat;  // Prospective Projection Matrix
  mat_4x4 rot_x;     // Rotation About X-Axis Matrix
  mat_4x4 rot_y;     // Rotation About Y-Axis Matrix
  mat_4x4 rot_z;     // Rotation About X-Axis Matrix
  mat_4x4 ident;     // Identity Matrix
  mat_4x4 cam_mat;   // Camera Effect Matrix

  // Fields
  float field_of_view;
  float zf;
  float zn;
  int screen_width;
  int screen_height;
  camera cam;
  std::vector<vec_3d> light_srcs;

  // Helper Math Functions
  void mat4x4_vec_mul(const vec_3d& vo, const mat_4x4& proj_mat, vec_3d& vp) {
    vp.x = vo.x * proj_mat.mat[0][0] + vo.y * proj_mat.mat[0][1] +
           vo.z * proj_mat.mat[0][2] + proj_mat.mat[0][3];
    vp.y = vo.x * proj_mat.mat[1][0] + vo.y * proj_mat.mat[1][1] +
           vo.z * proj_mat.mat[1][2] + proj_mat.mat[1][3];
    vp.z = vo.x * proj_mat.mat[2][0] + vo.y * proj_mat.mat[2][1] +
           vo.z * proj_mat.mat[2][2] + proj_mat.mat[2][3];
    vp.w = vo.x * proj_mat.mat[3][0] + vo.y * proj_mat.mat[3][1] +
           vo.z * proj_mat.mat[3][2] + proj_mat.mat[3][3];
  }

  float get_rad(float theta) { return (theta / 180.0f) * (22.0f / 7.0f); }

  void xy_scale(float scale_f, _3d_surface& os) {
    for (int i = 0; i < os.tris.size(); i++) {
      os.tris[i].p1.x *= (scale_f * this->screen_width);
      os.tris[i].p1.y *= (scale_f * this->screen_height);
      os.tris[i].p2.x *= (scale_f * this->screen_width);
      os.tris[i].p2.y *= (scale_f * this->screen_height);
      os.tris[i].p3.x *= (scale_f * this->screen_width);
      os.tris[i].p3.y *= (scale_f * this->screen_height);
    }
  }

  void pros_z_scale(_3d_surface& os) {
    for (int i = 0; i < os.tris.size(); i++) {
      this->scale_vec(os.tris[i].p1, (1.0f / os.tris[i].p1.w));
      this->scale_vec(os.tris[i].p2, (1.0f / os.tris[i].p2.w));
      this->scale_vec(os.tris[i].p3, (1.0f / os.tris[i].p3.w));
    }
  }

  void vec_minus(const vec_3d& p1, const vec_3d& p2, vec_3d& diff) {
    diff.x = p1.x - p2.x;
    diff.y = p1.y - p2.y;
    diff.z = p1.z - p2.z;
  }

  void vec_minus(const vec_2d& p1, const vec_2d& p2, vec_2d& diff) {
    diff.x = p1.x - p2.x;
    diff.y = p1.y - p2.y;
  }

  float get_vec_mag(const vec_3d& v1) {
    return sqrtf(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z);
  }

  float get_vec_mag(const vec_2d& v1) {
    return sqrtf(v1.x * v1.x + v1.y * v1.y);
  }

  void normalize_vec(vec_3d& v1) {
    float vec_mag = this->get_vec_mag(v1);
    this->scale_vec(v1, (1.0f / vec_mag));
  }

  void get_vec_norm(const vec_2d& v1, vec_2d& out_v) {
    float vec_mag = this->get_vec_mag(v1);
    out_v.x = v1.x / vec_mag;
    out_v.y = v1.y / vec_mag;
  }

  void get_norm(const triangle& surface, vec_3d& norm) {
    vec_3d v1;
    this->vec_minus(surface.p2, surface.p1, v1);
    vec_3d v2;
    this->vec_minus(surface.p3, surface.p1, v2);
    this->cross_prod(v1, v2, norm);
    this->normalize_vec(norm);
  }

  float dot_prod(const vec_3d& v1, const vec_3d& v2) {
    return (v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
  }

  void cross_prod(const vec_3d& v1, const vec_3d& v2, vec_3d& _out) {
    _out.x = v1.y * v2.z - v1.z * v2.y;
    _out.y = v1.z * v2.x - v1.x * v2.z;
    _out.z = v1.x * v2.y - v1.y * v2.x;
  }

  void scale_vec(vec_3d& v, float scale_f) {
    v.x *= scale_f;
    v.y *= scale_f;
    v.z *= scale_f;
  }

  void construct_cam_mat() {
    // Gram-Schmit Process to Get New Camera X Vector
    this->cam.x_vec.x -= (this->dot_prod(this->cam.x_vec, this->cam.z_vec) /
                          this->dot_prod(this->cam.z_vec, this->cam.z_vec)) *
                         this->cam.z_vec.x;
    this->cam.x_vec.y -= (this->dot_prod(this->cam.x_vec, this->cam.z_vec) /
                          this->dot_prod(this->cam.z_vec, this->cam.z_vec)) *
                         this->cam.z_vec.y;
    this->cam.x_vec.z -= (this->dot_prod(this->cam.x_vec, this->cam.z_vec) /
                          this->dot_prod(this->cam.z_vec, this->cam.z_vec)) *
                         this->cam.z_vec.z;
    this->cross_prod(this->cam.z_vec, this->cam.x_vec, this->cam.y_vec);

    this->cam_mat.mat[0][0] = this->cam.x_vec.x;
    this->cam_mat.mat[0][1] = this->cam.x_vec.y;
    this->cam_mat.mat[0][2] = this->cam.x_vec.z;
    this->cam_mat.mat[0][3] =
        -1.0f * this->dot_prod(this->cam.cam_vec, this->cam.x_vec);

    this->cam_mat.mat[1][0] = this->cam.y_vec.x;
    this->cam_mat.mat[1][1] = this->cam.y_vec.y;
    this->cam_mat.mat[1][2] = this->cam.y_vec.z;
    this->cam_mat.mat[1][3] =
        -1.0f * this->dot_prod(this->cam.cam_vec, this->cam.y_vec);

    this->cam_mat.mat[2][0] = this->cam.z_vec.x;
    this->cam_mat.mat[2][1] = this->cam.z_vec.y;
    this->cam_mat.mat[2][2] = this->cam.z_vec.z;
    this->cam_mat.mat[2][3] =
        -1.0f * this->dot_prod(this->cam.cam_vec, this->cam.z_vec);

    this->cam_mat.mat[3][3] = 1;
  }

  // Test if Certain Point is Front The Plane Face or Not
  bool is_pfp(const vec_3d& plane_p0,
              const vec_3d& plane_norm,
              const vec_3d& target_point) {
    vec_3d v;
    v.x = target_point.x - plane_p0.x;
    v.y = target_point.y - plane_p0.y;
    v.z = target_point.z - plane_p0.z;
    return (this->dot_prod(v, plane_norm) > 0.0f);
  }

  // Get Intersection of Plane With Line
  float get_plane_line_p(const vec_3d& plane_p0,
                         const vec_3d& plane_norm,
                         const vec_3d& line_p1,
                         const vec_3d& line_p2,
                         vec_3d& out_p) {
    vec_3d u;
    u.x = line_p2.x - line_p1.x;
    u.y = line_p2.y - line_p1.y;
    u.z = line_p2.z - line_p1.z;
    float c = (this->dot_prod(plane_p0, plane_norm) -
               this->dot_prod(line_p1, plane_norm)) /
              this->dot_prod(u, plane_norm);
    out_p.x = line_p1.x + c * u.x;
    out_p.y = line_p1.y + c * u.y;
    out_p.z = line_p1.z + c * u.z;

    vec_3d inner_vec;
    vec_3d entire_vec;
    this->vec_minus(out_p, line_p1, inner_vec);
    this->vec_minus(line_p2, line_p1, entire_vec);
    float inner_vec_mag = this->get_vec_mag(inner_vec);
    float entire_vec_mag = this->get_vec_mag(entire_vec);
    float inter_ratio = inner_vec_mag / entire_vec_mag;
    return inter_ratio;
  }

 public:
  zdx_core(float view_angel,
           float world_depth,
           int _screen_width,
           int _screen_height,
           camera _cam_data,
           std::vector<vec_3d> _light_srcs,
           float _zn)
      : field_of_view((view_angel / 180.0f) * (22.0f / 7.0f)),
        zf(world_depth),
        screen_width(_screen_width),
        screen_height(_screen_height),
        zn(_zn),
        cam(_cam_data),
        light_srcs(_light_srcs) {
    float aspect_ratio = (float)this->screen_height / this->screen_width;
    float plane_scale = 1 / tanf(this->field_of_view / 2);
    float z_scale = zf / (zf - zn);
// Debug Information
#ifdef DEBUG_MODE_PARAMS
    cout << "Aspect Ratio -> " << aspect_ratio << endl;
    cout << "Screen Width -> " << this->screen_height << endl;
    cout << "Screen Height -> " << this->screen_width << endl;
    cout << "CAM -> { " << this->camv.x << ", " << this->camv.y << ", "
         << this->camv.z << " }" << endl;
#endif
    // Constructing Transformation Matrices
    this->proj_mat.mat[0][0] = aspect_ratio * plane_scale;
    this->proj_mat.mat[1][1] = plane_scale;
    this->proj_mat.mat[2][2] = z_scale;
    this->proj_mat.mat[3][2] = 1.0;
    this->proj_mat.mat[2][3] = -zn * z_scale;

    this->rot_x.mat[0][0] = 1;
    this->rot_x.mat[3][3] = 1;

    this->rot_y.mat[1][1] = 1;
    this->rot_y.mat[3][3] = 1;

    this->rot_z.mat[2][2] = 1;
    this->rot_z.mat[3][3] = 1;

    this->ident.mat[0][0] = 1;
    this->ident.mat[1][1] = 1;
    this->ident.mat[2][2] = 1;
    this->ident.mat[3][3] = 1;

    // Camera Construction
    this->cross_prod(this->cam.z_vec, this->cam.x_vec, this->cam.y_vec);
    this->construct_cam_mat();
  }

  void x_axis_rot(float theta, const _3d_surface& os, _3d_surface& ts) {
    float theta_rad = this->get_rad(theta);
    this->rot_x.mat[1][1] = cosf(theta_rad);
    this->rot_x.mat[1][2] = -sinf(theta_rad);
    this->rot_x.mat[2][1] = sinf(theta_rad);
    this->rot_x.mat[2][2] = cosf(theta_rad);
    for (int i = 0; i < os.tris.size(); i++) {
      mat4x4_vec_mul(os.tris[i].p1, this->rot_x, ts.tris[i].p1);
      mat4x4_vec_mul(os.tris[i].p2, this->rot_x, ts.tris[i].p2);
      mat4x4_vec_mul(os.tris[i].p3, this->rot_x, ts.tris[i].p3);
    }
  }

  void y_axis_rot(float theta, const _3d_surface& os, _3d_surface& ts) {
    float theta_rad = this->get_rad(theta);
    this->rot_y.mat[0][0] = cosf(theta_rad);
    this->rot_y.mat[0][2] = -sinf(theta_rad);
    this->rot_y.mat[2][0] = sinf(theta_rad);
    this->rot_y.mat[2][2] = cosf(theta_rad);
    for (int i = 0; i < os.tris.size(); i++) {
      mat4x4_vec_mul(os.tris[i].p1, this->rot_y, ts.tris[i].p1);
      mat4x4_vec_mul(os.tris[i].p2, this->rot_y, ts.tris[i].p2);
      mat4x4_vec_mul(os.tris[i].p3, this->rot_y, ts.tris[i].p3);
    }
  }

  void z_axis_rot(float theta, const _3d_surface& os, _3d_surface& ts) {
    float theta_rad = this->get_rad(theta);
    this->rot_z.mat[0][0] = cosf(theta_rad);
    this->rot_z.mat[0][1] = -sinf(theta_rad);
    this->rot_z.mat[1][0] = sinf(theta_rad);
    this->rot_z.mat[1][1] = cosf(theta_rad);
    for (int i = 0; i < os.tris.size(); i++) {
      mat4x4_vec_mul(os.tris[i].p1, this->rot_z, ts.tris[i].p1);
      mat4x4_vec_mul(os.tris[i].p2, this->rot_z, ts.tris[i].p2);
      mat4x4_vec_mul(os.tris[i].p3, this->rot_z, ts.tris[i].p3);
    }
  }

  void x_translate(float shift, _3d_surface& os) {
    for (int i = 0; i < os.tris.size(); i++) {
      os.tris[i].p1.x += shift;
      os.tris[i].p2.x += shift;
      os.tris[i].p3.x += shift;
    }
  }

  void y_translate(float shift, _3d_surface& os) {
    for (int i = 0; i < os.tris.size(); i++) {
      os.tris[i].p1.y += shift;
      os.tris[i].p2.y += shift;
      os.tris[i].p3.y += shift;
    }
  }

  void z_translate(float shift, _3d_surface& os) {
    for (int i = 0; i < os.tris.size(); i++) {
      os.tris[i].p1.z += shift;
      os.tris[i].p2.z += shift;
      os.tris[i].p3.z += shift;
    }
  }

  void project(const _3d_surface& os, _3d_surface& ts, float scale_f) {
    for (int i = 0; i < os.tris.size(); i++) {
      mat4x4_vec_mul(os.tris[i].p1, this->proj_mat, ts.tris[i].p1);
      mat4x4_vec_mul(os.tris[i].p2, this->proj_mat, ts.tris[i].p2);
      mat4x4_vec_mul(os.tris[i].p3, this->proj_mat, ts.tris[i].p3);
    }
    // Scale All Vertices By Z Value To Generate Prospective Illusion
    this->pros_z_scale(ts);
    // Make Objects Bigger So We Can See, Because Projection Matrix Produce Very
    // Small Values
    this->xy_scale(scale_f, ts);
    // Move The World Origin To The Center of The Screen
    this->x_translate(this->screen_width / 2, ts);
    this->y_translate(this->screen_height / 2, ts);
  }

  void cam_view_filter(const _3d_surface& os, _3d_surface& fs) {
    this->construct_cam_mat();
    for (int i = 0; i < os.tris.size(); i++) {
      vec_3d surface_norm;
      this->get_norm(os.tris[i], surface_norm);
      vec_3d cam_ray;
      this->vec_minus(this->cam.cam_vec, os.tris[i].p1, cam_ray);
      if (this->dot_prod(surface_norm, cam_ray) > 0.0f) {
        triangle temp = os.tris[i];
        this->mat4x4_vec_mul(os.tris[i].p1, this->cam_mat, temp.p1);
        this->mat4x4_vec_mul(os.tris[i].p2, this->cam_mat, temp.p2);
        this->mat4x4_vec_mul(os.tris[i].p3, this->cam_mat, temp.p3);
        fs.tris.push_back(temp);
      }
    }
  }

  void illuminate(_3d_surface& os) {
    for (int i = 0; i < os.tris.size(); i++) {
      vec_3d surface_norm;
      this->get_norm(os.tris[i], surface_norm);
      // Apply One Direction Illumination
      vec_3d ill_vec = this->light_srcs[0];
      this->normalize_vec(ill_vec);
      this->scale_vec(ill_vec, -1.0f);
      float dot_prod = this->dot_prod(surface_norm, ill_vec);
      float gray_shade = dot_prod * 255;
      if (!(os.tris[i].face_color == olc::RED ||
            os.tris[i].face_color == olc::GREEN ||
            os.tris[i].face_color == olc::BLUE)) {
        os.tris[i].face_color = olc::Pixel(gray_shade, gray_shade, gray_shade);
      }
    }
  }

  void z_index(_3d_surface& os) {
    sort(os.tris.begin(), os.tris.end(), [](triangle& t1, triangle& t2) {
      float z1 = (t1.p1.z + t1.p2.z + t1.p3.z) / 3.0f;
      float z2 = (t2.p1.z + t2.p2.z + t2.p3.z) / 3.0f;
      return z1 > z2;
    });
  }

  void mat4x4_mat4x4_mul(const mat_4x4& mat1,
                         const mat_4x4& mat2,
                         mat_4x4& out_mat) {
    for (int i = 0; i < 4; i++) {
      vec_3d temp;
      mat4x4_vec_mul(
          {mat2.mat[0][i], mat2.mat[1][i], mat2.mat[2][i], mat2.mat[3][i]},
          mat1, temp);
      out_mat.mat[0][i] = temp.x;
      out_mat.mat[1][i] = temp.y;
      out_mat.mat[2][i] = temp.z;
      out_mat.mat[3][i] = temp.w;
    }
  }

  const mat_4x4& get_rotx_mat(float theta) {
    float theta_rad = this->get_rad(theta);
    this->rot_x.mat[1][1] = cosf(theta_rad);
    this->rot_x.mat[1][2] = -sinf(theta_rad);
    this->rot_x.mat[2][1] = sinf(theta_rad);
    this->rot_x.mat[2][2] = cosf(theta_rad);
    return this->rot_x;
  }

  const mat_4x4& get_roty_mat(float theta) {
    float theta_rad = this->get_rad(theta);
    this->rot_y.mat[0][0] = cosf(theta_rad);
    this->rot_y.mat[0][2] = -sinf(theta_rad);
    this->rot_y.mat[2][0] = sinf(theta_rad);
    this->rot_y.mat[2][2] = cosf(theta_rad);
    return this->rot_y;
  }

  const mat_4x4& get_rotz_mat(float theta) {
    float theta_rad = this->get_rad(theta);
    this->rot_z.mat[0][0] = cosf(theta_rad);
    this->rot_z.mat[0][1] = -sinf(theta_rad);
    this->rot_z.mat[1][0] = sinf(theta_rad);
    this->rot_z.mat[1][1] = cosf(theta_rad);
    return this->rot_z;
  }

  void general_transform(const mat_4x4& trans_mat,
                         const _3d_surface& os,
                         _3d_surface& ts) {
    for (int i = 0; i < os.tris.size(); i++) {
      mat4x4_vec_mul(os.tris[i].p1, trans_mat, ts.tris[i].p1);
      mat4x4_vec_mul(os.tris[i].p2, trans_mat, ts.tris[i].p2);
      mat4x4_vec_mul(os.tris[i].p3, trans_mat, ts.tris[i].p3);
    }
  }

  const mat_4x4& get_id_mat() { return this->ident; }

  camera& get_camera() { return this->cam; }

  void yaw_cam(float theta) {
    const mat_4x4& yaw_mat = this->get_roty_mat(theta);
    vec_3d temp;
    this->mat4x4_vec_mul(this->cam.z_vec, yaw_mat, temp);
    this->cam.z_vec.x = temp.x;
    this->cam.z_vec.y = temp.y;
    this->cam.z_vec.z = temp.z;
    this->cam.z_vec.w = temp.w;
  }

  void world_clip(const _3d_surface& os,
                  _3d_surface& _out_s,
                  const vec_3d& plane_point,
                  const vec_3d& plane_norm) {
    for (int i = 0; i < os.tris.size(); i++) {
      int pfpn = 0;
      std::vector<vec_3d> in_points, out_points;
      std::vector<vec_2d> in_tpoints, out_tpoints;
      if (this->is_pfp(plane_point, plane_norm, os.tris[i].p1)) {
        pfpn += 1;
        in_points.push_back(os.tris[i].p1);
        in_tpoints.push_back(os.tris[i].tp1);
      } else {
        out_points.push_back(os.tris[i].p1);
        out_tpoints.push_back(os.tris[i].tp1);
      }
      if (this->is_pfp(plane_point, plane_norm, os.tris[i].p2)) {
        pfpn += 1;
        in_points.push_back(os.tris[i].p2);
        in_tpoints.push_back(os.tris[i].tp2);
      } else {
        out_points.push_back(os.tris[i].p2);
        out_tpoints.push_back(os.tris[i].tp2);
      }
      if (this->is_pfp(plane_point, plane_norm, os.tris[i].p3)) {
        pfpn += 1;
        in_points.push_back(os.tris[i].p3);
        in_tpoints.push_back(os.tris[i].tp3);
      } else {
        out_points.push_back(os.tris[i].p3);
        out_tpoints.push_back(os.tris[i].tp3);
      }

      if (pfpn == 3) {
        _out_s.tris.push_back(os.tris[i]);
      } else if (pfpn == 2) {
        vec_3d int_p1;
        vec_2d int_tp1;
        float inter_ratio1 = this->get_plane_line_p(
            plane_point, plane_norm, in_points[0], out_points[0], int_p1);
        vec_3d int_p2;
        vec_2d int_tp2;
        float inter_ratio2 = this->get_plane_line_p(
            plane_point, plane_norm, in_points[1], out_points[0], int_p2);
        if (os.has_texture) {
          vec_2d tex_vec1;
          vec_2d tex_vec2;
          this->vec_minus(out_tpoints[0], in_tpoints[0], tex_vec1);
          this->vec_minus(out_tpoints[0], in_tpoints[1], tex_vec2);
          int_tp1.x = in_tpoints[0].x + inter_ratio1 * tex_vec1.x;
          int_tp1.y = in_tpoints[0].y + inter_ratio1 * tex_vec1.y;
          int_tp2.x = in_tpoints[1].x + inter_ratio2 * tex_vec2.x;
          int_tp2.y = in_tpoints[1].y + inter_ratio2 * tex_vec2.y;
        }
        triangle ct1, ct2;
        ct1.p1 = {in_points[0].x, in_points[0].y, in_points[0].z};
        ct1.p2 = {int_p1.x, int_p1.y, int_p1.z};
        ct1.p3 = {int_p2.x, int_p2.y, int_p2.z};
        ct1.tp1 = {in_tpoints[0].x, in_tpoints[0].y};
        ct1.tp2 = {int_tp1.x, int_tp1.y};
        ct1.tp3 = {int_tp2.x, int_tp2.y};
#ifdef DEBUG_MODE_GFX
        ct1.face_color = olc::Pixel(255, 0, 0);
#else
        ct1.face_color = os.tris[i].face_color;
#endif
        ct2.p1 = {in_points[0].x, in_points[0].y, in_points[0].z};
        ct2.p2 = {int_p2.x, int_p2.y, int_p2.z};
        ct2.p3 = {in_points[1].x, in_points[1].y, in_points[1].z};
        ct2.tp1 = {in_tpoints[0].x, in_tpoints[0].y};
        ct2.tp2 = {int_tp2.x, int_tp2.y};
        ct2.tp3 = {in_tpoints[1].x, in_tpoints[1].y};
#ifdef DEBUG_MODE_GFX
        ct2.face_color = olc::Pixel(0, 255, 0);
#else
        ct2.face_color = os.tris[i].face_color;
#endif
        _out_s.tris.push_back(ct1);
        _out_s.tris.push_back(ct2);
      } else if (pfpn == 1) {
        vec_3d int_p1;
        vec_2d int_tp1;
        float inter_ratio1 = this->get_plane_line_p(
            plane_point, plane_norm, in_points[0], out_points[0], int_p1);
        vec_3d int_p2;
        vec_2d int_tp2;
        float inter_ratio2 = this->get_plane_line_p(
            plane_point, plane_norm, in_points[0], out_points[1], int_p2);
        if (os.has_texture) {
          vec_2d tex_vec1;
          vec_2d tex_vec2;
          this->vec_minus(out_tpoints[0], in_tpoints[0], tex_vec1);
          this->vec_minus(out_tpoints[1], in_tpoints[0], tex_vec2);
          int_tp1.x = in_tpoints[0].x + inter_ratio1 * tex_vec1.x;
          int_tp1.y = in_tpoints[0].y + inter_ratio1 * tex_vec1.y;
          int_tp2.x = in_tpoints[0].x + inter_ratio2 * tex_vec2.x;
          int_tp2.y = in_tpoints[0].y + inter_ratio2 * tex_vec2.y;
        }
        triangle ct;
        ct.p1 = {in_points[0].x, in_points[0].y, in_points[0].z};
        ct.p2 = {int_p1.x, int_p1.y, int_p1.z};
        ct.p3 = {int_p2.x, int_p2.y, int_p2.z};
        ct.tp1 = {in_tpoints[0].x, in_tpoints[0].y};
        ct.tp2 = {int_tp1.x, int_tp1.y};
        ct.tp3 = {int_tp2.x, int_tp2.y};
#ifdef DEBUG_MODE_GFX
        ct.face_color = olc::Pixel(0, 0, 255);
#else
        ct.face_color = os.tris[i].face_color;
#endif
        _out_s.tris.push_back(ct);
      }
    }
  }

  void screen_clip(const _3d_surface& os, _3d_surface& _out_s) {
    vec_3d planes_points[4], planes_norms[4];

    planes_points[0] = {0.0f, this->screen_height / 2.0f, 0.0f};
    planes_points[1] = {this->screen_width / 2.0f, (float)this->screen_height,
                        0.0f};
    planes_points[2] = {(float)this->screen_width, this->screen_height / 2.0f,
                        0.0f};
    planes_points[3] = {this->screen_width / 2.0f, 0.0f, 0.0f};

    planes_norms[0] = {1.0f, 0.0f, 0.0f};
    planes_norms[1] = {0.0f, -1.0f, 0.0f};
    planes_norms[2] = {-1.0f, 0.0f, 0.0f};
    planes_norms[3] = {0.0f, 1.0f, 0.0f};

    for (int i = 0; i < os.tris.size(); i++) {
      std::vector<triangle> tris_to_clip;
      tris_to_clip.push_back(os.tris[i]);
      for (int j = 0; j < 4; j++) {
        std::vector<triangle> new_tris;
        for (int k = 0; k < tris_to_clip.size(); k++) {
          int pfpn = 0;
          triangle tri_to_clip = tris_to_clip[k];
          std::vector<vec_3d> in_points, out_points;
          std::vector<vec_2d> in_tpoints, out_tpoints;
          if (this->is_pfp(planes_points[j], planes_norms[j], tri_to_clip.p1)) {
            pfpn += 1;
            in_points.push_back(tri_to_clip.p1);
            in_tpoints.push_back(tri_to_clip.tp1);
          } else {
            out_points.push_back(tri_to_clip.p1);
            out_tpoints.push_back(tri_to_clip.tp1);
          }
          if (this->is_pfp(planes_points[j], planes_norms[j], tri_to_clip.p2)) {
            pfpn += 1;
            in_points.push_back(tri_to_clip.p2);
            in_tpoints.push_back(tri_to_clip.tp2);
          } else {
            out_points.push_back(tri_to_clip.p2);
            out_tpoints.push_back(tri_to_clip.tp2);
          }
          if (this->is_pfp(planes_points[j], planes_norms[j], tri_to_clip.p3)) {
            pfpn += 1;
            in_points.push_back(tri_to_clip.p3);
            in_tpoints.push_back(tri_to_clip.tp3);
          } else {
            out_points.push_back(tri_to_clip.p3);
            out_tpoints.push_back(tri_to_clip.tp3);
          }

          if (pfpn == 3) {
            new_tris.push_back(tri_to_clip);
          } else if (pfpn == 2) {
            vec_3d int_p1;
            vec_2d int_tp1;
            float inter_ratio1 =
                this->get_plane_line_p(planes_points[j], planes_norms[j],
                                       in_points[0], out_points[0], int_p1);
            vec_3d int_p2;
            vec_2d int_tp2;
            float inter_ratio2 =
                this->get_plane_line_p(planes_points[j], planes_norms[j],
                                       in_points[1], out_points[0], int_p2);
            if (os.has_texture) {
              vec_2d tex_vec1;
              vec_2d tex_vec2;
              this->vec_minus(out_tpoints[0], in_tpoints[0], tex_vec1);
              this->vec_minus(out_tpoints[0], in_tpoints[1], tex_vec2);
              int_tp1.x = in_tpoints[0].x + inter_ratio1 * tex_vec1.x;
              int_tp1.y = in_tpoints[0].y + inter_ratio1 * tex_vec1.y;
              int_tp2.x = in_tpoints[1].x + inter_ratio2 * tex_vec2.x;
              int_tp2.y = in_tpoints[1].y + inter_ratio2 * tex_vec2.y;
            }
            triangle ct1, ct2;
            ct1.p1 = {in_points[0].x, in_points[0].y, in_points[0].z};
            ct1.p2 = {int_p1.x, int_p1.y, int_p1.z};
            ct1.p3 = {int_p2.x, int_p2.y, int_p2.z};
            ct1.tp1 = {in_tpoints[0].x, in_tpoints[0].y};
            ct1.tp2 = {int_tp1.x, int_tp1.y};
            ct1.tp3 = {int_tp2.x, int_tp2.y};
#ifdef DEBUG_MODE_GFX
            ct1.face_color = olc::Pixel(255, 0, 0);
#else
            ct1.face_color = os.tris[i].face_color;
#endif
            ct2.p1 = {in_points[0].x, in_points[0].y, in_points[0].z};
            ct2.p2 = {int_p2.x, int_p2.y, int_p2.z};
            ct2.p3 = {in_points[1].x, in_points[1].y, in_points[1].z};
            ct2.tp1 = {in_tpoints[0].x, in_tpoints[0].y};
            ct2.tp2 = {int_tp2.x, int_tp2.y};
            ct2.tp3 = {in_tpoints[1].x, in_tpoints[1].y};
#ifdef DEBUG_MODE_GFX
            ct2.face_color = olc::Pixel(0, 255, 0);
#else
            ct2.face_color = os.tris[i].face_color;
#endif
            new_tris.push_back(ct1);
            new_tris.push_back(ct2);
          } else if (pfpn == 1) {
            vec_3d int_p1;
            vec_2d int_tp1;
            float inter_ratio1 =
                this->get_plane_line_p(planes_points[j], planes_norms[j],
                                       in_points[0], out_points[0], int_p1);
            vec_3d int_p2;
            vec_2d int_tp2;
            float inter_ratio2 =
                this->get_plane_line_p(planes_points[j], planes_norms[j],
                                       in_points[0], out_points[1], int_p2);
            if (os.has_texture) {
              vec_2d tex_vec1;
              vec_2d tex_vec2;
              this->vec_minus(out_tpoints[0], in_tpoints[0], tex_vec1);
              this->vec_minus(out_tpoints[1], in_tpoints[0], tex_vec2);
              int_tp1.x = in_tpoints[0].x + inter_ratio1 * tex_vec1.x;
              int_tp1.y = in_tpoints[0].y + inter_ratio1 * tex_vec1.y;
              int_tp2.x = in_tpoints[0].x + inter_ratio2 * tex_vec2.x;
              int_tp2.y = in_tpoints[0].y + inter_ratio2 * tex_vec2.y;
            }
            triangle ct;
            ct.p1 = {in_points[0].x, in_points[0].y, in_points[0].z};
            ct.p2 = {int_p1.x, int_p1.y, int_p1.z};
            ct.p3 = {int_p2.x, int_p2.y, int_p2.z};
            ct.tp1 = {in_tpoints[0].x, in_tpoints[0].y};
            ct.tp2 = {int_tp1.x, int_tp1.y};
            ct.tp3 = {int_tp2.x, int_tp2.y};
#ifdef DEBUG_MODE_GFX
            ct.face_color = olc::Pixel(0, 0, 255);
#else
            ct.face_color = os.tris[i].face_color;
#endif
            new_tris.push_back(ct);
          }
        }
        tris_to_clip.clear();
        for (int k = 0; k < new_tris.size(); k++) {
          triangle temp;
          temp.p1 = {new_tris[k].p1.x, new_tris[k].p1.y, new_tris[k].p1.z};
          temp.p2 = {new_tris[k].p2.x, new_tris[k].p2.y, new_tris[k].p2.z};
          temp.p3 = {new_tris[k].p3.x, new_tris[k].p3.y, new_tris[k].p3.z};
          temp.tp1 = {new_tris[k].tp1.x, new_tris[k].tp1.y};
          temp.tp2 = {new_tris[k].tp2.x, new_tris[k].tp2.y};
          temp.tp3 = {new_tris[k].tp3.x, new_tris[k].tp3.y};
          temp.face_color = new_tris[k].face_color;
          tris_to_clip.push_back(temp);
        }
      }
      for (int j = 0; j < tris_to_clip.size(); j++) {
        _out_s.tris.push_back(tris_to_clip[j]);
      }
    }
  }

  _3d_surface bulk_pipeline(const _3d_surface& os, float _theta = 0.0f) {
    // 3D Transformation
    _3d_surface ts = os;  // 3D Transformed Surface
    // mat_4x4 _3d_trans_mat;
    // this->mat4x4_mat4x4_mul(this->get_rotz_mat(_theta),
    //                         this->get_rotx_mat(_theta), _3d_trans_mat);
    // this->general_transform(_3d_trans_mat, os, ts);
    this->z_translate(4.0f, ts);
    // Filter Only Front Traingles to Apply The Transformation on
    _3d_surface front_ts;
    front_ts.has_texture = ts.has_texture;
    this->cam_view_filter(ts, front_ts);
    // Clipping Against Plane Front of The Camera
    _3d_surface world_clipped;
    world_clipped.has_texture = front_ts.has_texture;
    this->world_clip(front_ts, world_clipped, {0.0f, 0.0f, 2.0f},
                     {0.0f, 0.0f, 1.0f});  // Clip Plane is Front of the Camera
    // Adding Illumination
    // this->illuminate(world_clipped);
    // Projection
    _3d_surface front_ts_proj = world_clipped;
    this->project(world_clipped, front_ts_proj, 0.5f);
    _3d_surface front_ts_proj_screen_clipped;
    this->screen_clip(front_ts_proj, front_ts_proj_screen_clipped);
    this->z_index(front_ts_proj_screen_clipped);
    return front_ts_proj_screen_clipped;
  }
};

class zdx_texture {
 private:
  olc::Pixel* raw_data;
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
      this->raw_data = new olc::Pixel[this->tex_width * this->tex_height];
      for (int y = 0; y < this->tex_height; y++) {
        png_bytep row = row_pointers[y];
        for (int x = 0; x < this->tex_width; x++) {
          png_bytep px = &(row[x * 4]);
          this->raw_data[y * this->tex_width + x] =
              olc::Pixel(px[0], px[1], px[2]);
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

  olc::Pixel get_pixel(int x, int y) const {
    if (x >= 0 && x < this->tex_width && y >= 0 && y < this->tex_height) {
      return this->raw_data[y * this->tex_width + x];
    } else {
      return olc::Pixel(255, 0, 0);
    }
  }

  olc::Pixel get_pixel_norm(float norm_x, float norm_y) const {
    return this->get_pixel(norm_x * this->tex_width, norm_y * this->tex_height);
  }
};
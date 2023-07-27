// Copyright 2021 Adobe Research. All rights reserved.
// To view a copy of the license, visit LICENSE.md.

// Copyright (C) 2008 Gael Guennebaud <gael.guennebaud@inria.fr>
//
// Eigen is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// Eigen is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License or the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// Eigen. If not, see <http://www.gnu.org/licenses/>.

#ifndef CAMERA_H
#define CAMERA_H

#if EIGEN_COMP_MSVC
// NOTE MSVC often gives C4127 warnings with compiletime if statements. See bug
// 1362. This workaround is ugly, but it does the job.
#define EIGEN_CONST_CONDITIONAL(cond) (void)0, cond
#else
#define EIGEN_CONST_CONDITIONAL(cond) cond
#endif

#include "common.h"
#include "io/serialization.h"

#include <map>

/// Represents a 3D frame, i.e., an orthogonal basis with the position of the
/// origin.
class Frame {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  inline Frame(const Vector3f &pos = Vector3f::Zero(),
               const Quaternionf &o = Quaternionf())
      : orientation(o), position(pos) {}
  Frame lerp(real_t alpha, const Frame &other) const {
    return Frame((1.f - alpha) * position + alpha * other.position,
                 orientation.slerp(alpha, other.orientation));
  }

  Matrix4f getMatrix() const;

  Quaternionf orientation;
  Vector3f position;
};

class KeyFrame {
public:
  KeyFrame() {}
  KeyFrame(const Quaternionf &q, const Vector3f &pos);

  Vector3f position() const { return p_; }
  Quaternionf orientation() const { return q_; }
  Vector3f tgP() const { return tgP_; }
  Quaternionf tgQ() const { return tgQ_; }
  void computeTangent(const KeyFrame *const prev, const KeyFrame *const next);

  static Frame lerpFrame(real_t alpha, const KeyFrame *a, const KeyFrame *b);
  static Quaternionf squad(const Quaternionf &a,
                                  const Quaternionf &tgA,
                                  const Quaternionf &tgB,
                                  const Quaternionf &b, real_t t);

private:
  Vector3f p_, tgP_;
  Quaternionf q_, tgQ_;
};

class Camera;

class KeyFrameInterpolator {
public:

  KeyFrameInterpolator()
      : mAlpha(0), mFramerate(30.f), mInterpLinear(true),
        mTangentAreValid(false), mLoop(true) {}

  void setCamera(Camera *camera);
  void addKeyFrame(int t, const Quaternionf &q,
                   const Vector3f &po);
  void addKeyFrame(int delta_t = 30);
  Frame interpolatedFrame();
  void reset(const Frame &initFrame);
  void setInterpolationMode(bool on);

  real_t getTime() const { return mAlpha; }
  void setTime(real_t t);

  int gotoPrevKeyFrame();
  int gotoNextKeyFrame();
  void deleteKeyFrame();

  int update(real_t delta_t);

  using TimeLine = std::map<int, std::unique_ptr<KeyFrame>>;
  const TimeLine &timeline() const { return mTimeline; }

private:
  TimeLine mTimeline;
  real_t mAlpha;
  real_t mFramerate;
  bool mInterpLinear;
  bool mTangentAreValid;
  bool mLoop;
  Camera *mCamera;

};

/// Represents a virtual camera, which is essentially a Frame with a given view
/// frustum and viewport
class Camera {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  Camera();
  Camera(const Camera &other);
  Camera &operator=(const Camera &other);

  //------------------------------------------------------------
  // Viewport setters and getters
  void setViewport(unsigned int offsetx, unsigned int offsety,
                   unsigned int width, unsigned int height);
  void setViewport(unsigned int width, unsigned int height);
  inline unsigned int vpX() const { return mVpX; }
  inline unsigned int vpY() const { return mVpY; }
  inline unsigned int vpWidth() const { return mVpWidth; }
  inline unsigned int vpHeight() const { return mVpHeight; }
  //------------------------------------------------------------

  //------------------------------------------------------------
  // View frustum setters and getters
  /// \returns the vertical field of view angle (in radian)
  inline real_t fovY() const { return mFovY; }
  /// sets the vertical field of view angle (in radian)
  void setFovY(real_t value);
  /// \returns the distance of the image (near) plane
  real_t nearDist() const { return mNearDist; }
  real_t farDist() const { return mFarDist; }
  /** Setup the perspective projection matrix */
  void setPerspective(real_t fovY, real_t near, real_t far);
  //------------------------------------------------------------

  //------------------------------------------------------------
  // Frame setters and getters
  /// sets the position of the camera
  void setPosition(const Vector3f &pos);
  /// \returns the position of the camera
  inline const Vector3f &position() const { return mFrame.position; }
  /// sets the orientation of the camera
  void setOrientation(const Quaternionf &q);
  /// \returns the orientation of the camera
  inline const Quaternionf &orientation() const {
    return mFrame.orientation;
  }
  void setFrame(const Frame &f);
  const Frame &frame(void) const { return mFrame; }
  /// sets the view direction of the camera
  void setDirection(const Vector3f &newDirection);
  /// \returns the view direction, i.e., the -z axis of the frame
  Vector3f direction() const;
  /// \returns the view distance
  real_t dist() const;
  /// \returns the up vertical direction, i.e., the y axis of the frame
  Vector3f up() const;
  /// \returns the right horizontal direction , i.e., the x axis of the frame
  Vector3f right() const;
  //------------------------------------------------------------

  //------------------------------------------------------------
  // Advanced Frame setters
  /** Setup the camera position and orientation based on its \a position, \a a
   * target point, \a the up vector */
  void lookAt(const Vector3f &position, const Vector3f &target,
              const Vector3f &up);
  /// \returns the priviligied view target point
  inline const Vector3f &target(void) const { return mTarget; }
  /// sets the target of the camera
  void setTarget(const Vector3f &target);
  //------------------------------------------------------------

  /// \returns the affine transformation from camera to global space
  const Affine3f &viewMatrix() const;
  /// \returns the projective transformation matrix from camera space to the
  /// normalized image space
  const Matrix4f &projectionMatrix() const;

  /// rotates the camera around the target point using the rotation \a q
  void rotateAroundTarget(const Quaternionf &q);
  /// rotates the camera around the own camera position using the rotation \a q
  void localRotate(const Quaternionf &q);
  /// moves the camera toward the target
  void zoom(real_t d);
  /// moves the camera toward the target
  void zoomByFactor(real_t f);
  /// moves the camera by \a t defined in the local camera space
  void localTranslate(const Vector3f &t);

  Vector3f unProject(const Vector2f &uv, real_t depth,
                            const Matrix4f &invModelview) const;
  /// project a given point from the image space to the global space
  Vector3f unProject(const Vector2f &uv, real_t depth) const;

  KeyFrameInterpolator &getInterpolator() { return mKeyFrameInterpolators; }

  void setInterpolatedFrame(Frame f);

  bool load(const std::string &filename, int& startframe, int& endframe);
  void save(const std::string &filename) const;

  const Matrix3Xf &getPoints() const { return mPoints; }

protected:
  void updateViewMatrix() const;
  void updateProjectionMatrix() const;

  unsigned int mVpX, mVpY;
  unsigned int mVpWidth, mVpHeight;

  Frame mFrame;
  int mStartframe, mEndframe;

  mutable Affine3f mViewMatrix;
  mutable Matrix4f mProjectionMatrix;

  mutable bool mViewIsUptodate;
  mutable bool mProjIsUptodate;

  // used by rotateAroundTarget
  Vector3f mTarget;

  real_t mFovY;
  real_t mNearDist;
  real_t mFarDist;

  mutable Matrix3Xf mPoints;

  KeyFrameInterpolator mKeyFrameInterpolators;

  contess_befriend_serialization(Camera);
};

#endif // CAMERA_H

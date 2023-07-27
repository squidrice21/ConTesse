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

#include "camera.h"
#include <fstream>
#include <iostream>

Matrix4f Frame::getMatrix() const {
  Affine3f transformation;
  Quaternionf q = orientation.conjugate();
  transformation.linear() = q.toRotationMatrix();
  transformation.translation() = -(transformation.linear() * position);
  Matrix4f M = transformation.matrix();
  return M.inverse();
}

KeyFrame::KeyFrame(const Quaternionf &q, const Vector3f &pos) {
  p_ = pos;
  q_ = q;
}

// generic linear interpolation method
template <typename T> T lerp(real_t t, const T &a, const T &b) {
  return a * (1 - t) + b * t;
}

// quaternion slerp
Quaternionf slerp(real_t t, const Quaternionf &a, const Quaternionf &b,
                  bool allowFlip = true) {
  static const real_t one = real_t(1) - Eigen::NumTraits<real_t>::epsilon();
  real_t d = a.dot(b);
  real_t absD = abs(d);

  real_t scale0;
  real_t scale1;

  if (absD >= one) {
    scale0 = real_t(1) - t;
    scale1 = t;
  } else {
    real_t theta = acos(absD);
    real_t sinTheta = sin(theta);
    scale0 = sin((real_t(1) - t) * theta) / sinTheta;
    scale1 = sin((t * theta)) / sinTheta;
  }
  if (allowFlip && d < 0)
    scale1 = -scale1;

  return Quaternionf(scale0 * a.coeffs() + scale1 * b.coeffs());
}

/*! Returns the slerp interpolation of the two Quaternions \p a and \p b, at
  time \p t, using tangents \p tgA and \p tgB. The resulting Quaternion is
  "between" \p a and \p b (result is \p a when \p t=0 and \p b for \p t=1).
*/
Quaternionf KeyFrame::squad(const Quaternionf &a, const Quaternionf &tgA,
                            const Quaternionf &tgB, const Quaternionf &b,
                            real_t t) {
  Quaternionf ab = slerp(t, a, b, true);
  Quaternionf tg = slerp(t, tgA, tgB, false);
  return slerp(2.0 * t * (1.0 - t), ab, tg, false);
}

Frame KeyFrame::lerpFrame(real_t alpha, const KeyFrame *a, const KeyFrame *b) {
  return Frame(lerp(alpha, a->position(), b->position()),
               Quaternionf(slerp(alpha, a->orientation(), b->orientation())));
}

/*! Returns the logarithm of the Quaternion. See also exp(). */
inline Quaternionf log(const Quaternionf &q) {
  real_t len = sqrt(q.x() * q.x() + q.y() * q.y() + q.z() * q.z());

  if (len < 1E-6)
    return Quaternionf(0.0, q.x(), q.y(), q.z());
  else {
    real_t coef = acos(q.w()) / len;
    return Quaternionf(0.0, q.x() * coef, q.y() * coef, q.z() * coef);
  }
}

/*! Returns the exponential of the Quaternion. See also log(). */
inline Quaternionf exp(const Quaternionf &q) {
  real_t theta = sqrt(q.x() * q.x() + q.y() * q.y() + q.z() * q.z());

  if (theta < 1E-6)
    return Quaternionf(cos(theta), q.x(), q.y(), q.z());
  else {
    real_t coef = sin(theta) / theta;
    return Quaternionf(cos(theta), q.x() * coef, q.y() * coef, q.z() * coef);
  }
}

/*! Returns log(a. inverse() * b). Useful for squadTangent(). */
inline Quaternionf lnDif(const Quaternionf &a, const Quaternionf &b) {
  Quaternionf dif = a.inverse() * b;
  dif.normalize();
  return ::log(dif);
}

/*! Returns a tangent Quaternion for \p center, defined by \p before and \p
 after Quaternions. Useful for smooth spline interpolation of Quaternion with
 squad() and slerp(). */
inline Quaternionf squadTangent(const Quaternionf &before,
                                const Quaternionf &center,
                                const Quaternionf &after) {
  Quaternionf l1 = ::lnDif(center, before);
  Quaternionf l2 = ::lnDif(center, after);
  Quaternionf e = Quaternionf(-0.25 * (l1.coeffs() + l2.coeffs()));
  e = center * ::exp(e);

  return e;
}

void KeyFrame::computeTangent(const KeyFrame *const prev,
                              const KeyFrame *const next) {
  tgP_ = 0.5 * (next->position() - prev->position());
  tgQ_ = ::squadTangent(prev->orientation(), q_, next->orientation());
}

void KeyFrameInterpolator::setCamera(Camera *camera) { mCamera = camera; }

void KeyFrameInterpolator::addKeyFrame(int t, const Quaternionf &q,
                                       const Vector3f &po) {
  mTimeline[t] = std::make_unique<KeyFrame>(q, po);
}

void KeyFrameInterpolator::addKeyFrame(int delta_t) {
  int t;
  if (mTimeline.empty()) {
    t = 0;
  } else {
    t = (--mTimeline.end())->first + delta_t;
  }
  mTimeline[t] =
      std::make_unique<KeyFrame>(Quaternionf(mCamera->viewMatrix().linear()),
                                 mCamera->viewMatrix().translation());
  mTangentAreValid = false;
}

int KeyFrameInterpolator::update(real_t delta_t) {
  if (mTimeline.empty())
    return 0;
  mAlpha += delta_t;
  mCamera->setInterpolatedFrame(interpolatedFrame());
  return int(mAlpha * mFramerate);
}

void KeyFrameInterpolator::setTime(real_t t) {
  mAlpha = t / mFramerate;
  mCamera->setInterpolatedFrame(interpolatedFrame());
}

int KeyFrameInterpolator::gotoPrevKeyFrame() {
  TimeLine::const_iterator lo = mTimeline.upper_bound(int(mAlpha * mFramerate));
  if (lo != mTimeline.begin())
    --lo;
  if (lo->first == int(mAlpha * mFramerate) && lo != mTimeline.begin())
    --lo;
  setTime(real_t(lo->first));
  return lo->first;
}

int KeyFrameInterpolator::gotoNextKeyFrame() {
  TimeLine::const_iterator hi = mTimeline.upper_bound(int(mAlpha * mFramerate));
  setTime(real_t(hi->first));
  return hi->first;
}

void KeyFrameInterpolator::deleteKeyFrame() {
  TimeLine::const_iterator lo = mTimeline.lower_bound(int(mAlpha * mFramerate));
  mTimeline.erase(lo->first);
}

Frame KeyFrameInterpolator::interpolatedFrame() {
  TimeLine::const_iterator hi = mTimeline.upper_bound(int(mAlpha * mFramerate));
  TimeLine::const_iterator lo = hi;
  if (lo != mTimeline.begin())
    --lo;

  Frame currentFrame;

  if (hi == mTimeline.end()) {
    if (mLoop) {
      mAlpha = 0.f;
      currentFrame = Frame(mTimeline.begin()->second->position(),
                           mTimeline.begin()->second->orientation());
    } else {
      // end
      currentFrame = Frame(lo->second->position(), lo->second->orientation());
    }
  } else if (hi == mTimeline.begin()) {
    // start
    currentFrame = Frame(hi->second->position(), hi->second->orientation());
  } else {
    real_t s = (mAlpha * mFramerate - lo->first) / (hi->first - lo->first);

    if (mInterpLinear) {
      currentFrame = KeyFrame::lerpFrame(s, lo->second.get(), hi->second.get());
    } else {
      if (!mTangentAreValid) {
        mTangentAreValid = true;
        // Update tangents
        TimeLine::iterator prev = mTimeline.begin();
        TimeLine::iterator kf = prev;
        while (kf != mTimeline.end()) {
          TimeLine::iterator next = prev;
          next++;
          if (next != mTimeline.end())
            kf->second->computeTangent(prev->second.get(), next->second.get());
          else
            kf->second->computeTangent(prev->second.get(), kf->second.get());
          prev = kf;
          kf = next;
        }
      }

      Vector3f delta = hi->second->position() - lo->second->position();
      Vector3f v1 = 3.0 * delta - 2.0 * lo->second->tgP() - hi->second->tgP();
      Vector3f v2 = -2.0 * delta + lo->second->tgP() + hi->second->tgP();

      Vector3f pos =
          lo->second->position() + s * (lo->second->tgP() + s * (v1 + s * v2));
      Quaternionf q =
          KeyFrame::squad(lo->second->orientation(), lo->second->tgQ(),
                          hi->second->tgQ(), hi->second->orientation(), s);
      currentFrame = Frame(pos, q);
    }

    currentFrame.orientation.coeffs().normalize();
  }

  currentFrame.orientation = currentFrame.orientation.inverse();
  currentFrame.position = -(currentFrame.orientation * currentFrame.position);

  return currentFrame;
}

void KeyFrameInterpolator::reset(const Frame &initFrame) {
  mTimeline.clear();

  Frame aux0 = mCamera->frame();
  aux0.orientation = aux0.orientation.inverse();
  aux0.position = mCamera->viewMatrix().translation();
  mTimeline[0] = std::make_unique<KeyFrame>(aux0.orientation, aux0.position);

  mCamera->setTarget(Vector3f::Zero());

  // compute the rotation duration to move the camera to the target
  Frame aux1 = mCamera->frame();
  aux1.orientation = aux1.orientation.inverse();
  aux1.position = mCamera->viewMatrix().translation();
  real_t duration = aux0.orientation.angularDistance(aux1.orientation) * 0.9;
  if (duration < 0.1)
    duration = 0.1;

  // put the camera at that time step:
  aux1 = aux0.lerp(duration / 2, initFrame);
  // and make it look at the target again
  aux1.orientation = aux1.orientation.inverse();
  aux1.position = -(aux1.orientation * aux1.position);
  mCamera->setFrame(aux1);
  mCamera->setTarget(Vector3f::Zero());

  // add this camera keyframe
  aux1.orientation = aux1.orientation.inverse();
  aux1.position = mCamera->viewMatrix().translation();
  mTimeline[duration] =
      std::make_unique<KeyFrame>(aux1.orientation, aux1.position);

  mTangentAreValid = false;

  mAlpha = 0;
  interpolatedFrame();
}

void KeyFrameInterpolator::setInterpolationMode(bool on) { mInterpLinear = on; }

Camera::Camera() : mViewIsUptodate(false), mProjIsUptodate(false) {
  mViewMatrix.setIdentity();

  mFovY = M_PI / 3.;
  mNearDist = 0.1;
  mFarDist = 50000.;

  mVpX = 0;
  mVpY = 0;

  setPosition(Vector3f::Constant(100.));
  setTarget(Vector3f::Zero());
  mKeyFrameInterpolators.setCamera(this);
}

Camera &Camera::operator=(const Camera &other) {
  mVpX = other.mVpX;
  mVpY = other.mVpY;
  mVpWidth = other.mVpWidth;
  mVpHeight = other.mVpHeight;

  mFrame = other.mFrame;

  mViewMatrix = other.viewMatrix();
  mProjectionMatrix = other.projectionMatrix();

  mViewIsUptodate = true;
  mProjIsUptodate = true;

  mTarget = other.mTarget;

  mFovY = other.mFovY;
  mNearDist = other.mNearDist;
  mFarDist = other.mFarDist;

  mPoints = other.mPoints;

  return *this;
}

Camera::Camera(const Camera &other) { *this = other; }

void Camera::setPerspective(real_t fovY, real_t near, real_t far) {
  mFovY = fovY;
  mNearDist = near;
  mFarDist = far;
  mProjIsUptodate = false;
}

void Camera::setViewport(unsigned int offsetx, unsigned int offsety,
                         unsigned int width, unsigned int height) {
  mVpX = offsetx;
  mVpY = offsety;
  mVpWidth = width;
  mVpHeight = height;

  mProjIsUptodate = false;
}

void Camera::setViewport(unsigned int width, unsigned int height) {
  mVpWidth = width;
  mVpHeight = height;

  mProjIsUptodate = false;
}

void Camera::setFovY(real_t value) {
  mFovY = value;
  mProjIsUptodate = false;
}

void Camera::setDirection(const Vector3f &newDirection) {
  Vector3f up = Vector3f(0, 0, 1); // this->up();

  Matrix3f camAxes;

  camAxes.col(2) = (-newDirection).normalized();
  camAxes.col(0) = up.cross(camAxes.col(2)).normalized();
  camAxes.col(1) = camAxes.col(2).cross(camAxes.col(0)).normalized();
  setOrientation(Quaternionf(camAxes));

  mViewIsUptodate = false;
}

Vector3f Camera::direction() const {
  updateViewMatrix();
  return -mViewMatrix.linear().row(2);
}

real_t Camera::dist() const { return (position() - mTarget).norm(); }
Vector3f Camera::up() const {
  updateViewMatrix();
  return mViewMatrix.linear().row(1);
}
Vector3f Camera::right() const {
  updateViewMatrix();
  return mViewMatrix.linear().row(0);
}

void Camera::lookAt(const Vector3f &position, const Vector3f &target,
                    const Vector3f &up) {
  mTarget = target;
  mFrame.position = position;
  Matrix3f R;
  R.col(2) = (position - target).normalized();
  R.col(0) = up.cross(R.col(2)).normalized();
  R.col(1) = R.col(2).cross(R.col(0));
  setOrientation(Quaternionf(R));
  mViewIsUptodate = false;
}

void Camera::setPosition(const Vector3f &p) {
  mFrame.position = p;
  mViewIsUptodate = false;
}

void Camera::setTarget(const Vector3f &target) {
  mTarget = target;
  if (!mTarget.isApprox(position())) {
    Vector3f newDirection = mTarget - position();
    setDirection(newDirection.normalized());
  }
}

void Camera::setOrientation(const Quaternionf &q) {
  mFrame.orientation = q;
  mViewIsUptodate = false;
}

void Camera::setFrame(const Frame &f) {
  mFrame = f;
  mViewIsUptodate = false;
}

void Camera::rotateAroundTarget(const Quaternionf &q) {
  Matrix4f mrot, mt, mtm;

  // update the transform matrix
  updateViewMatrix();
  Vector3f t = mViewMatrix * mTarget;

  mViewMatrix = Translation3f(t) * q * Translation3f(-t) * mViewMatrix;

  Quaternionf qa(mViewMatrix.linear());
  qa = qa.conjugate();
  setOrientation(qa);
  setPosition(-(qa * mViewMatrix.translation()));

  mViewIsUptodate = true;
}

void Camera::localRotate(const Quaternionf &q) {
  real_t dist = (position() - mTarget).norm();
  setOrientation(orientation() * q);
  mTarget = position() + dist * direction();
  mViewIsUptodate = false;
}

void Camera::zoom(real_t d) {
  real_t dist = (position() - mTarget).norm();
  if (dist > d) {
    setPosition(position() + direction() * d);
    mViewIsUptodate = false;
  }
}

void Camera::zoomByFactor(real_t f) {
  real_t d = (position() - mTarget).norm() * f;
  setPosition(position() + direction() * d);
  mViewIsUptodate = false;
}

void Camera::localTranslate(const Vector3f &t) {
  Vector3f trans = orientation() * t;
  setPosition(position() + trans);
  mTarget += trans;

  mViewIsUptodate = false;
}

void Camera::updateViewMatrix() const {
  if (!mViewIsUptodate) {
    Quaternionf q = orientation().conjugate();
    mViewMatrix.linear() = q.toRotationMatrix();
    mViewMatrix.translation() = -(mViewMatrix.linear() * position());

    mViewIsUptodate = true;
  }
}

const Affine3f &Camera::viewMatrix() const {
  updateViewMatrix();
  return mViewMatrix;
}

void Camera::updateProjectionMatrix() const {
  if (!mProjIsUptodate) {
    mProjectionMatrix.setIdentity();
    real_t aspect = real_t(mVpWidth) / real_t(mVpHeight);
    real_t theta = mFovY * 0.5;
    real_t range = mFarDist - mNearDist;
    real_t invtan = 1. / tan(theta);

    mProjectionMatrix(0, 0) = invtan / aspect;
    mProjectionMatrix(1, 1) = invtan;
    mProjectionMatrix(2, 2) = -(mNearDist + mFarDist) / range;
    mProjectionMatrix(3, 2) = -1;
    mProjectionMatrix(2, 3) = -2 * mNearDist * mFarDist / range;
    mProjectionMatrix(3, 3) = 0;

    // pyramide
    real_t ym = tan(mFovY * 0.5);
    real_t xm = ((real_t)mVpWidth) * (ym * 1.0 / mVpHeight);
    real_t zm = 0.75f;
    mPoints = Matrix3Xf(3, 28);
    mPoints.col(0) = Vector3f(0.0, 0.0, 0.0);
    mPoints.col(1) = Vector3f(xm, ym, -zm);
    mPoints.col(2) = Vector3f(xm, ym, -zm);
    mPoints.col(3) = Vector3f(xm, -ym, -zm);
    mPoints.col(4) = Vector3f(xm, -ym, -zm);
    mPoints.col(5) = Vector3f(0.0, 0.0, 0.0);
    mPoints.col(6) = Vector3f(0.0, 0.0, 0.0);
    mPoints.col(7) = Vector3f(-xm, ym, -zm);
    mPoints.col(8) = Vector3f(-xm, ym, -zm);
    mPoints.col(9) = Vector3f(-xm, -ym, -zm);
    mPoints.col(10) = Vector3f(-xm, -ym, -zm);
    mPoints.col(11) = Vector3f(0.0, 0.0, 0.0);
    mPoints.col(12) = Vector3f(0.0, 0.0, 0.0);
    mPoints.col(13) = Vector3f(xm, ym, -zm);
    mPoints.col(14) = Vector3f(xm, ym, -zm);
    mPoints.col(15) = Vector3f(-xm, ym, -zm);
    mPoints.col(16) = Vector3f(-xm, ym, -zm);
    mPoints.col(17) = Vector3f(0.0, 0.0, 0.0);
    mPoints.col(18) = Vector3f(0.0, 0.0, 0.0);
    mPoints.col(19) = Vector3f(xm, -ym, -zm);
    mPoints.col(20) = Vector3f(xm, -ym, -zm);
    mPoints.col(21) = Vector3f(-xm, -ym, -zm);
    mPoints.col(22) = Vector3f(-xm, -ym, -zm);
    mPoints.col(23) = Vector3f(0.0, 0.0, 0.0);
    mPoints.col(24) = Vector3f(-xm, ym, -zm);
    mPoints.col(25) = Vector3f(0.0, 1.5 * ym, -zm);
    mPoints.col(26) = Vector3f(xm, ym, -zm);
    mPoints.col(27) = Vector3f(0.0, 1.5 * ym, -zm);

    mProjIsUptodate = true;
  }
}

const Matrix4f &Camera::projectionMatrix() const {
  updateProjectionMatrix();
  return mProjectionMatrix;
}

Vector3f Camera::unProject(const Vector2f &uv, real_t depth) const {
  Matrix4f inv = mViewMatrix.inverse().matrix();
  return unProject(uv, depth, inv);
}

Vector3f Camera::unProject(const Vector2f &uv, real_t depth,
                           const Matrix4f &invModelview) const {
  updateViewMatrix();
  updateProjectionMatrix();

  Vector3f a(2. * uv.x() / real_t(mVpWidth) - 1.,
             2. * uv.y() / real_t(mVpHeight) - 1., 1.);
  a.x() *= depth / mProjectionMatrix(0, 0);
  a.y() *= depth / mProjectionMatrix(1, 1);
  a.z() = -depth;
  // Vector4f b = invModelview * Vector4f(a.x(), a.y(), a.z(), 1.);
  // Vector3f b = mViewMatrix.linear().transpose() * a +
  // position();//Vector4f(a.x(), a.y(), a.z(), 1.);

  Vector3f b =
      a.x() * right() + a.y() * up() - a.z() * direction() + position();

  return Vector3f(b.x(), b.y(), b.z());
}

void Camera::setInterpolatedFrame(Frame f) { setFrame(f); }

bool Camera::load(const std::string &filename, int &startframe, int &endframe) {
  using namespace std;
  ifstream file;
  file.open(filename);
  if (!file.is_open()) {
    std::cout << "ERROR: Unable to open camera file " << filename << std::endl;
    return false;
  }

  int frame;
  float focal, sensor_height, sensor_width;
  Eigen::Quaternionf q;
  Eigen::Vector3f p;
  string line, prefix;
  while (getline(file, line)) {
    prefix = line.substr(0, 2);
    if (prefix == "sf") {
      sscanf(line.c_str(), "sf %d", &startframe);
    } else if (prefix == "ef") {
      sscanf(line.c_str(), "ef %d", &endframe);
    } else if (prefix == "fo") {
      sscanf(line.c_str(), "fo %f", &focal);
    } else if (prefix == "sv") {
      sscanf(line.c_str(), "sv %f", &sensor_width);
    } else if (prefix == "sh") {
      sscanf(line.c_str(), "sh %f", &sensor_height);
    } else {
      // frame number
      sscanf(line.c_str(), "%d %f %f %f %f %f %f %f", &frame, &p.x(), &p.y(),
             &p.z(), &q.w(), &q.x(), &q.y(), &q.z());
      q = q.inverse();
      p = -(q * p);
      mKeyFrameInterpolators.addKeyFrame(frame, q.cast<real_t>(),
                                         p.cast<real_t>());
    }
  }

  setViewport(0, 0, sensor_width, sensor_height);
  mFovY = focal;
  mStartframe = startframe;
  mEndframe = endframe;

  file.close();
  return true;
}

void Camera::save(const std::string &filename) const {
  const KeyFrameInterpolator::TimeLine &timeline =
      mKeyFrameInterpolators.timeline();

  if (timeline.size() < 1) {
    return;
  }

  using namespace std;
  ofstream file;
  file.open(filename, ofstream::out);
  if (!file.is_open()) {
    std::cout << "ERROR: Unable to open file " << filename << std::endl;
  }

  file << "sf " << (*timeline.begin()).first << endl;
  file << "ef " << (*(--timeline.end())).first << endl;
  file << "fo " << mFovY << endl;
  file << "sv " << mVpWidth << endl;
  file << "sh " << mVpHeight << endl;

  for (const auto &key : timeline) {
    Quaternionf q = key.second->orientation();
    Vector3f p = key.second->position();
    q = q.inverse();
    p = -(q * p);
    file << key.first << " " << p.x() << " " << p.y() << " " << p.z() << " "
         << q.w() << " " << q.x() << " " << q.y() << " " << q.z() << endl;
  }

  file.close();
}

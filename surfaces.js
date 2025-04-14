/**
 * surfaces.js
 * Definitions for the 2D functions (surfaces) used in gradient descent demos.
 */

const surfaces = {
  singleWell: {
    name: "Single Well (Convex Bowl)",
    eq: "$f(w_0, w_1) = w_0^2 + w_1^2$",
    range: [-3, 3],
    step: 0.1,
    f: (w0, w1) => (w0*w0 + w1*w1),
    grad: (w0, w1) => [2*w0, 2*w1]
  },
  doubleWell: {
    name: "Double Well",
    eq: "$f(w_0, w_1) = (w_0^2 - 1)^2 + w_1^2$",
    range: [-2, 2],
    step: 0.05,
    f: (w0, w1) => {
      const part = (w0*w0 - 1);
      return part*part + w1*w1;
    },
    grad: (w0, w1) => {
      const gw0 = 4*w0*(w0*w0 - 1); // derivative wrt w0
      const gw1 = 2*w1;            // derivative wrt w1
      return [gw0, gw1];
    }
  },
  fourWell: {
    name: "Four Well",
    eq: "$f(w_0, w_1) = (w_0^2 - 1)^2 + (w_1^2 - 1)^2$",
    range: [-2, 2],
    step: 0.05,
    f: (w0, w1) => {
      const p0 = (w0*w0 - 1);
      const p1 = (w1*w1 - 1);
      return p0*p0 + p1*p1;
    },
    grad: (w0, w1) => {
      const gw0 = 4*w0*(w0*w0 - 1);
      const gw1 = 4*w1*(w1*w1 - 1);
      return [gw0, gw1];
    }
  }
}; 
// kriging-contour v1.0.2 Copyright 2020 freegis
(function (global, factory) {
typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
typeof define === 'function' && define.amd ? define(['exports'], factory) :
(global = global || self, factory(global.kriging = global.kriging || {}));
}(this, (function (exports) { 'use strict';

function ascending(a, b) {
  return a < b ? -1 : a > b ? 1 : a >= b ? 0 : NaN;
}

function bisector(compare) {
  if (compare.length === 1) compare = ascendingComparator(compare);
  return {
    left: function(a, x, lo, hi) {
      if (lo == null) lo = 0;
      if (hi == null) hi = a.length;
      while (lo < hi) {
        var mid = lo + hi >>> 1;
        if (compare(a[mid], x) < 0) lo = mid + 1;
        else hi = mid;
      }
      return lo;
    },
    right: function(a, x, lo, hi) {
      if (lo == null) lo = 0;
      if (hi == null) hi = a.length;
      while (lo < hi) {
        var mid = lo + hi >>> 1;
        if (compare(a[mid], x) > 0) hi = mid;
        else lo = mid + 1;
      }
      return lo;
    }
  };
}

function ascendingComparator(f) {
  return function(d, x) {
    return ascending(f(d), x);
  };
}

var ascendingBisect = bisector(ascending);

function extent(values, valueof) {
  var n = values.length,
      i = -1,
      value,
      min,
      max;

  if (valueof == null) {
    while (++i < n) { // Find the first comparable value.
      if ((value = values[i]) != null && value >= value) {
        min = max = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = values[i]) != null) {
            if (min > value) min = value;
            if (max < value) max = value;
          }
        }
      }
    }
  }

  else {
    while (++i < n) { // Find the first comparable value.
      if ((value = valueof(values[i], i, values)) != null && value >= value) {
        min = max = value;
        while (++i < n) { // Compare the remaining values.
          if ((value = valueof(values[i], i, values)) != null) {
            if (min > value) min = value;
            if (max < value) max = value;
          }
        }
      }
    }
  }

  return [min, max];
}

function range(start, stop, step) {
  start = +start, stop = +stop, step = (n = arguments.length) < 2 ? (stop = start, start = 0, 1) : n < 3 ? 1 : +step;

  var i = -1,
      n = Math.max(0, Math.ceil((stop - start) / step)) | 0,
      range = new Array(n);

  while (++i < n) {
    range[i] = start + i * step;
  }

  return range;
}

var e10 = Math.sqrt(50),
    e5 = Math.sqrt(10),
    e2 = Math.sqrt(2);

function tickStep(start, stop, count) {
  var step0 = Math.abs(stop - start) / Math.max(0, count),
      step1 = Math.pow(10, Math.floor(Math.log(step0) / Math.LN10)),
      error = step0 / step1;
  if (error >= e10) step1 *= 10;
  else if (error >= e5) step1 *= 5;
  else if (error >= e2) step1 *= 2;
  return stop < start ? -step1 : step1;
}

function thresholdSturges(values) {
  return Math.ceil(Math.log(values.length) / Math.LN2) + 1;
}

var array = Array.prototype;

var slice = array.slice;

function ascending$1(a, b) {
  return a - b;
}

function area(ring) {
  var i = 0, n = ring.length, area = ring[n - 1][1] * ring[0][0] - ring[n - 1][0] * ring[0][1];
  while (++i < n) area += ring[i - 1][1] * ring[i][0] - ring[i - 1][0] * ring[i][1];
  return area;
}

function constant(x) {
  return function() {
    return x;
  };
}

function contains(ring, hole) {
  var i = -1, n = hole.length, c;
  while (++i < n) if (c = ringContains(ring, hole[i])) return c;
  return 0;
}

function ringContains(ring, point) {
  var x = point[0], y = point[1], contains = -1;
  for (var i = 0, n = ring.length, j = n - 1; i < n; j = i++) {
    var pi = ring[i], xi = pi[0], yi = pi[1], pj = ring[j], xj = pj[0], yj = pj[1];
    if (segmentContains(pi, pj, point)) return 0;
    if (((yi > y) !== (yj > y)) && ((x < (xj - xi) * (y - yi) / (yj - yi) + xi))) contains = -contains;
  }
  return contains;
}

function segmentContains(a, b, c) {
  var i; return collinear(a, b, c) && within(a[i = +(a[0] === b[0])], c[i], b[i]);
}

function collinear(a, b, c) {
  return (b[0] - a[0]) * (c[1] - a[1]) === (c[0] - a[0]) * (b[1] - a[1]);
}

function within(p, q, r) {
  return p <= q && q <= r || r <= q && q <= p;
}

function noop() {}

var cases = [
  [],
  [[[1.0, 1.5], [0.5, 1.0]]],
  [[[1.5, 1.0], [1.0, 1.5]]],
  [[[1.5, 1.0], [0.5, 1.0]]],
  [[[1.0, 0.5], [1.5, 1.0]]],
  [[[1.0, 1.5], [0.5, 1.0]], [[1.0, 0.5], [1.5, 1.0]]],
  [[[1.0, 0.5], [1.0, 1.5]]],
  [[[1.0, 0.5], [0.5, 1.0]]],
  [[[0.5, 1.0], [1.0, 0.5]]],
  [[[1.0, 1.5], [1.0, 0.5]]],
  [[[0.5, 1.0], [1.0, 0.5]], [[1.5, 1.0], [1.0, 1.5]]],
  [[[1.5, 1.0], [1.0, 0.5]]],
  [[[0.5, 1.0], [1.5, 1.0]]],
  [[[1.0, 1.5], [1.5, 1.0]]],
  [[[0.5, 1.0], [1.0, 1.5]]],
  []
];

function d3_contours() {
  var dx = 1,
      dy = 1,
      threshold = thresholdSturges,
      smooth = smoothLinear;

  function contours(values) {
    var tz = threshold(values);

    // Convert number of thresholds into uniform thresholds.
    if (!Array.isArray(tz)) {
      var domain = extent(values), start = domain[0], stop = domain[1];
      tz = tickStep(start, stop, tz);
      tz = range(Math.floor(start / tz) * tz, Math.floor(stop / tz) * tz, tz);
    } else {
      tz = tz.slice().sort(ascending$1);
    }

    return tz.map(function(value) {
      return contour(values, value);
    });
  }

  // Accumulate, smooth contour rings, assign holes to exterior rings.
  // Based on https://github.com/mbostock/shapefile/blob/v0.6.2/shp/polygon.js
  function contour(values, value) {
    var polygons = [],
        holes = [];

    isorings(values, value, function(ring) {
      smooth(ring, values, value);
      if (area(ring) > 0) polygons.push([ring]);
      else holes.push(ring);
    });

    holes.forEach(function(hole) {
      for (var i = 0, n = polygons.length, polygon; i < n; ++i) {
        if (contains((polygon = polygons[i])[0], hole) !== -1) {
          polygon.push(hole);
          return;
        }
      }
    });

    return {
      type: "MultiPolygon",
      value: value,
      coordinates: polygons
    };
  }

  // Marching squares with isolines stitched into rings.
  // Based on https://github.com/topojson/topojson-client/blob/v3.0.0/src/stitch.js
  function isorings(values, value, callback) {
    var fragmentByStart = new Array,
        fragmentByEnd = new Array,
        x, y, t0, t1, t2, t3;

    // Special case for the first row (y = -1, t2 = t3 = 0).
    x = y = -1;
    t1 = values[0] >= value;
    cases[t1 << 1].forEach(stitch);
    while (++x < dx - 1) {
      t0 = t1, t1 = values[x + 1] >= value;
      cases[t0 | t1 << 1].forEach(stitch);
    }
    cases[t1 << 0].forEach(stitch);

    // General case for the intermediate rows.
    while (++y < dy - 1) {
      x = -1;
      t1 = values[y * dx + dx] >= value;
      t2 = values[y * dx] >= value;
      cases[t1 << 1 | t2 << 2].forEach(stitch);
      while (++x < dx - 1) {
        t0 = t1, t1 = values[y * dx + dx + x + 1] >= value;
        t3 = t2, t2 = values[y * dx + x + 1] >= value;
        cases[t0 | t1 << 1 | t2 << 2 | t3 << 3].forEach(stitch);
      }
      cases[t1 | t2 << 3].forEach(stitch);
    }

    // Special case for the last row (y = dy - 1, t0 = t1 = 0).
    x = -1;
    t2 = values[y * dx] >= value;
    cases[t2 << 2].forEach(stitch);
    while (++x < dx - 1) {
      t3 = t2, t2 = values[y * dx + x + 1] >= value;
      cases[t2 << 2 | t3 << 3].forEach(stitch);
    }
    cases[t2 << 3].forEach(stitch);

    function stitch(line) {
      var start = [line[0][0] + x, line[0][1] + y],
          end = [line[1][0] + x, line[1][1] + y],
          startIndex = index(start),
          endIndex = index(end),
          f, g;
      if (f = fragmentByEnd[startIndex]) {
        if (g = fragmentByStart[endIndex]) {
          delete fragmentByEnd[f.end];
          delete fragmentByStart[g.start];
          if (f === g) {
            f.ring.push(end);
            callback(f.ring);
          } else {
            fragmentByStart[f.start] = fragmentByEnd[g.end] = {start: f.start, end: g.end, ring: f.ring.concat(g.ring)};
          }
        } else {
          delete fragmentByEnd[f.end];
          f.ring.push(end);
          fragmentByEnd[f.end = endIndex] = f;
        }
      } else if (f = fragmentByStart[endIndex]) {
        if (g = fragmentByEnd[startIndex]) {
          delete fragmentByStart[f.start];
          delete fragmentByEnd[g.end];
          if (f === g) {
            f.ring.push(end);
            callback(f.ring);
          } else {
            fragmentByStart[g.start] = fragmentByEnd[f.end] = {start: g.start, end: f.end, ring: g.ring.concat(f.ring)};
          }
        } else {
          delete fragmentByStart[f.start];
          f.ring.unshift(start);
          fragmentByStart[f.start = startIndex] = f;
        }
      } else {
        fragmentByStart[startIndex] = fragmentByEnd[endIndex] = {start: startIndex, end: endIndex, ring: [start, end]};
      }
    }
  }

  function index(point) {
    return point[0] * 2 + point[1] * (dx + 1) * 4;
  }

  function smoothLinear(ring, values, value) {
    ring.forEach(function(point) {
      var x = point[0],
          y = point[1],
          xt = x | 0,
          yt = y | 0,
          v0,
          v1 = values[yt * dx + xt];
      if (x > 0 && x < dx && xt === x) {
        v0 = values[yt * dx + xt - 1];
        point[0] = x + (value - v0) / (v1 - v0) - 0.5;
      }
      if (y > 0 && y < dy && yt === y) {
        v0 = values[(yt - 1) * dx + xt];
        point[1] = y + (value - v0) / (v1 - v0) - 0.5;
      }
    });
  }

  contours.contour = contour;

  contours.size = function(_) {
    if (!arguments.length) return [dx, dy];
    var _0 = Math.ceil(_[0]), _1 = Math.ceil(_[1]);
    if (!(_0 > 0) || !(_1 > 0)) throw new Error("invalid size");
    return dx = _0, dy = _1, contours;
  };

  contours.thresholds = function(_) {
    return arguments.length ? (threshold = typeof _ === "function" ? _ : Array.isArray(_) ? constant(slice.call(_)) : constant(_), contours) : threshold;
  };

  contours.smooth = function(_) {
    return arguments.length ? (smooth = _ ? smoothLinear : noop, contours) : smooth === smoothLinear;
  };

  return contours;
}

//数组最大值
Array.prototype.max = function () {
	return Math.max.apply(null, this);
};
//数组最小值
Array.prototype.min = function () {
	return Math.min.apply(null, this);
};
//数组平均值
Array.prototype.mean = function () {
	var i,
	sum;
	for (i = 0, sum = 0; i < this.length; i++)
		sum += this[i];
	return sum / this.length;
};

//将数组第一项取出为v，生成长度为n的数组，每个数组item为v
Array.prototype.rep = function (n) {
	var arrayn = new Array(n);
	var v = this[0];
	for (var i = 0; i < n; i++) {
		arrayn[i] = v;
	}
	return arrayn;
};

Array.prototype.pip = function (x, y) {
	var i,
	j,
	c = false;
	for (i = 0, j = this.length - 1; i < this.length; j = i++) {
		if (((this[i][1] > y) != (this[j][1] > y)) &&
			(x < (this[j][0] - this[i][0]) * (y - this[i][1]) / (this[j][1] - this[i][1]) + this[i][0])) {
			c = !c;
		}
	}
	return c;
};

// Matrix algebra
function kriging_matrix_diag(c, n) {
	var i,
	Z = [0].rep(n * n);
	for (i = 0; i < n; i++)
		Z[i * n + i] = c;
	return Z;
}function kriging_matrix_transpose(X, n, m) {
	var i,
	j,
	Z = Array(m * n);
	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
			Z[j * n + i] = X[i * m + j];
	return Z;
}function kriging_matrix_add(X, Y, n, m) {
	var i,
	j,
	Z = Array(n * m);
	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
			Z[i * m + j] = X[i * m + j] + Y[i * m + j];
	return Z;
}// Naive matrix multiplication
function kriging_matrix_multiply(X, Y, n, m, p) {
	var i,
	j,
	k,
	Z = Array(n * p);
	for (i = 0; i < n; i++) {
		for (j = 0; j < p; j++) {
			Z[i * p + j] = 0;
			for (k = 0; k < m; k++)
				Z[i * p + j] += X[i * m + k] * Y[k * p + j];
		}
	}
	return Z;
}// Cholesky decomposition
function kriging_matrix_chol(X, n) {
	var i,
	j,
	k,
	p = Array(n);
	for (i = 0; i < n; i++)
		p[i] = X[i * n + i];
	for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++)
			p[i] -= X[i * n + j] * X[i * n + j];
		if (p[i] <= 0)
			return false;
		p[i] = Math.sqrt(p[i]);
		for (j = i + 1; j < n; j++) {
			for (k = 0; k < i; k++)
				X[j * n + i] -= X[j * n + k] * X[i * n + k];
			X[j * n + i] /= p[i];
		}
	}
	for (i = 0; i < n; i++)
		X[i * n + i] = p[i];
	return true;
}// Inversion of cholesky decomposition
function kriging_matrix_chol2inv(X, n) {
	var i,
	j,
	k,
	sum;
	for (i = 0; i < n; i++) {
		X[i * n + i] = 1 / X[i * n + i];
		for (j = i + 1; j < n; j++) {
			sum = 0;
			for (k = i; k < j; k++)
				sum -= X[j * n + k] * X[k * n + i];
			X[j * n + i] = sum / X[j * n + j];
		}
	}
	for (i = 0; i < n; i++)
		for (j = i + 1; j < n; j++)
			X[i * n + j] = 0;
	for (i = 0; i < n; i++) {
		X[i * n + i] *= X[i * n + i];
		for (k = i + 1; k < n; k++)
			X[i * n + i] += X[k * n + i] * X[k * n + i];
		for (j = i + 1; j < n; j++)
			for (k = j; k < n; k++)
				X[i * n + j] += X[k * n + i] * X[k * n + j];
	}
	for (i = 0; i < n; i++)
		for (j = 0; j < i; j++)
			X[i * n + j] = X[j * n + i];

}// Inversion via gauss-jordan elimination
function kriging_matrix_solve(X, n) {
	var m = n;
	var b = Array(n * n);
	var indxc = Array(n);
	var indxr = Array(n);
	var ipiv = Array(n);
	var i,
	icol,
	irow,
	j,
	k,
	l,
	ll;
	var big,
	dum,
	pivinv,
	temp;

	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++) {
			if (i == j)
				b[i * n + j] = 1;
			else
				b[i * n + j] = 0;
		}
	for (j = 0; j < n; j++)
		ipiv[j] = 0;
	for (i = 0; i < n; i++) {
		big = 0;
		for (j = 0; j < n; j++) {
			if (ipiv[j] != 1) {
				for (k = 0; k < n; k++) {
					if (ipiv[k] == 0) {
						if (Math.abs(X[j * n + k]) >= big) {
							big = Math.abs(X[j * n + k]);
							irow = j;
							icol = k;
						}
					}
				}
			}
		}
		++(ipiv[icol]);

		if (irow != icol) {
			for (l = 0; l < n; l++) {
				temp = X[irow * n + l];
				X[irow * n + l] = X[icol * n + l];
				X[icol * n + l] = temp;
			}
			for (l = 0; l < m; l++) {
				temp = b[irow * n + l];
				b[irow * n + l] = b[icol * n + l];
				b[icol * n + l] = temp;
			}
		}
		indxr[i] = irow;
		indxc[i] = icol;

		if (X[icol * n + icol] == 0)
			return false; // Singular

		pivinv = 1 / X[icol * n + icol];
		X[icol * n + icol] = 1;
		for (l = 0; l < n; l++)
			X[icol * n + l] *= pivinv;
		for (l = 0; l < m; l++)
			b[icol * n + l] *= pivinv;

		for (ll = 0; ll < n; ll++) {
			if (ll != icol) {
				dum = X[ll * n + icol];
				X[ll * n + icol] = 0;
				for (l = 0; l < n; l++)
					X[ll * n + l] -= X[icol * n + l] * dum;
				for (l = 0; l < m; l++)
					b[ll * n + l] -= b[icol * n + l] * dum;
			}
		}
	}
	for (l = (n - 1); l >= 0; l--)
		if (indxr[l] != indxc[l]) {
			for (k = 0; k < n; k++) {
				temp = X[k * n + indxr[l]];
				X[k * n + indxr[l]] = X[k * n + indxc[l]];
				X[k * n + indxc[l]] = temp;
			}
		}

	return true;
}

// Variogram models
function kriging_variogram_gaussian(h, nugget, range, sill, A) {
	return nugget + ((sill - nugget) / range) *
	(1.0 - Math.exp( - (1.0 / A) * Math.pow(h / range, 2)));
}function kriging_variogram_exponential(h, nugget, range, sill, A) {
	return nugget + ((sill - nugget) / range) *
	(1.0 - Math.exp( - (1.0 / A) * (h / range)));
}function kriging_variogram_spherical(h, nugget, range, sill, A) {
	if (h > range)
		return nugget + (sill - nugget) / range;
	return nugget + ((sill - nugget) / range) *
	(1.5 * (h / range) - 0.5 * Math.pow(h / range, 3));
}
var kriging = {

};




// Train using gaussian processes with bayesian priors
kriging.train = function (t, x, y, model, sigma2, alpha) {
	var variogram = {
		t : t,
		x : x,
		y : y,
		nugget : 0.0,
		range : 0.0,
		sill : 0.0,
		A : 1 / 3,
		n : 0
	};
	switch (model) {
	case "gaussian":
		variogram.model = kriging_variogram_gaussian;
		break;
	case "exponential":
		variogram.model = kriging_variogram_exponential;
		break;
	case "spherical":
		variogram.model = kriging_variogram_spherical;
		break;
	}
	// Lag distance/semivariance
	var i,
	j,
	k,
	l,
	n = t.length;
	var distance = Array((n * n - n) / 2);
	for (i = 0, k = 0; i < n; i++)
		for (j = 0; j < i; j++, k++) {
			distance[k] = Array(2);
			distance[k][0] = Math.pow(
					Math.pow(x[i] - x[j], 2) +
					Math.pow(y[i] - y[j], 2), 0.5);
			distance[k][1] = Math.abs(t[i] - t[j]);
		}
	distance.sort(function (a, b) {
		return a[0] - b[0];
	});
	variogram.range = distance[(n * n - n) / 2 - 1][0];

	// Bin lag distance
	var lags = ((n * n - n) / 2) > 30 ? 30 : (n * n - n) / 2;
	var tolerance = variogram.range / lags;
	var lag = [0].rep(lags);
	var semi = [0].rep(lags);
	if (lags < 30) {
		for (l = 0; l < lags; l++) {
			lag[l] = distance[l][0];
			semi[l] = distance[l][1];
		}
	} else {
		for (i = 0, j = 0, k = 0, l = 0; i < lags && j < ((n * n - n) / 2); i++, k = 0) {
			while (distance[j][0] <= ((i + 1) * tolerance)) {
				lag[l] += distance[j][0];
				semi[l] += distance[j][1];
				j++;
				k++;
				if (j >= ((n * n - n) / 2))
					break;
			}
			if (k > 0) {
				lag[l] /= k;
				semi[l] /= k;
				l++;
			}
		}
		if (l < 2)
			return variogram; // Error: Not enough points
	}

	// Feature transformation
	n = l;
	variogram.range = lag[n - 1] - lag[0];
	var X = [1].rep(2 * n);
	var Y = Array(n);
	var A = variogram.A;
	for (i = 0; i < n; i++) {
		switch (model) {
		case "gaussian":
			X[i * 2 + 1] = 1.0 - Math.exp( - (1.0 / A) * Math.pow(lag[i] / variogram.range, 2));
			break;
		case "exponential":
			X[i * 2 + 1] = 1.0 - Math.exp( - (1.0 / A) * lag[i] / variogram.range);
			break;
		case "spherical":
			X[i * 2 + 1] = 1.5 * (lag[i] / variogram.range) -
				0.5 * Math.pow(lag[i] / variogram.range, 3);
			break;
		}		Y[i] = semi[i];
	}

	// Least squares
	var Xt = kriging_matrix_transpose(X, n, 2);
	var Z = kriging_matrix_multiply(Xt, X, 2, n, 2);
	Z = kriging_matrix_add(Z, kriging_matrix_diag(1 / alpha, 2), 2, 2);
	var cloneZ = Z.slice(0);
	if (kriging_matrix_chol(Z, 2))
		kriging_matrix_chol2inv(Z, 2);
	else {
		kriging_matrix_solve(cloneZ, 2);
		Z = cloneZ;
	}
	var W = kriging_matrix_multiply(kriging_matrix_multiply(Z, Xt, 2, 2, n), Y, 2, n, 1);

	// Variogram parameters
	variogram.nugget = W[0];
	variogram.sill = W[1] * variogram.range + variogram.nugget;
	variogram.n = x.length;

	// Gram matrix with prior
	n = x.length;
	var K = Array(n * n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) {
			K[i * n + j] = variogram.model(Math.pow(Math.pow(x[i] - x[j], 2) +
						Math.pow(y[i] - y[j], 2), 0.5),
					variogram.nugget,
					variogram.range,
					variogram.sill,
					variogram.A);
			K[j * n + i] = K[i * n + j];
		}
		K[i * n + i] = variogram.model(0, variogram.nugget,
				variogram.range,
				variogram.sill,
				variogram.A);
	}

	// Inverse penalized Gram matrix projected to target vector
	var C = kriging_matrix_add(K, kriging_matrix_diag(sigma2, n), n, n);
	var cloneC = C.slice(0);
	if (kriging_matrix_chol(C, n))
		kriging_matrix_chol2inv(C, n);
	else {
		kriging_matrix_solve(cloneC, n);
		C = cloneC;
	}

	// Copy unprojected inverted matrix as K
	var K = C.slice(0);
	var M = kriging_matrix_multiply(C, t, n, n, 1);
	variogram.K = K;
	variogram.M = M;

	return variogram;
};

// Model prediction
kriging.predict = function (x, y, variogram) {
	var i,
	k = Array(variogram.n);
	for (i = 0; i < variogram.n; i++)
		k[i] = variogram.model(Math.pow(Math.pow(x - variogram.x[i], 2) +
					Math.pow(y - variogram.y[i], 2), 0.5),
				variogram.nugget, variogram.range,
				variogram.sill, variogram.A);
	return kriging_matrix_multiply(k, variogram.M, 1, variogram.n, 1)[0];
};
kriging.variance = function (x, y, variogram) {
	var i,
	k = Array(variogram.n);
	for (i = 0; i < variogram.n; i++)
		k[i] = variogram.model(Math.pow(Math.pow(x - variogram.x[i], 2) +
					Math.pow(y - variogram.y[i], 2), 0.5),
				variogram.nugget, variogram.range,
				variogram.sill, variogram.A);
	return variogram.model(0, variogram.nugget, variogram.range,
		variogram.sill, variogram.A) +
	kriging_matrix_multiply(kriging_matrix_multiply(k, variogram.K,
			1, variogram.n, variogram.n),
		k, 1, variogram.n, 1)[0];
};

// 生成克里金grid
kriging.getGridInfo = function (bbox,variogram,width) {
	var grid = [];
	//x方向
	var xlim=[bbox[0],bbox[2]];
	var ylim=[bbox[1],bbox[3]];
	var zlim=[variogram.t.min(), variogram.t.max()];
	
	//xy方向地理跨度
	var geoX_width=xlim[1]-xlim[0];
	var geoY_width=ylim[1]-ylim[0];

	//如果x_width设置，初始基于200计算。
	let x_width,y_width;
	if(!width)
		x_width=200;
	else
		x_width=Math.ceil(width);
	//让图像的xy比例与地理的xy比例保持一致
	y_width=Math.ceil(x_width/(geoX_width/geoY_width));

	//地理跨度/图像跨度=当前地图图上分辨率
	var x_resolution=geoX_width*1.0/x_width;
	var y_resolution=geoY_width*1.0/y_width;
	
	var xtarget,ytarget;

	for (let j = 0; j < y_width; j++) {
		for (let k =0; k <x_width; k++) {
			xtarget = bbox[0] + k * x_resolution;
			ytarget = bbox[1] + j * y_resolution;
			grid.push(kriging.predict(xtarget, ytarget, variogram));
		}
	}
	return {
		grid : grid,
		n : x_width,
		m : y_width,
		xlim : xlim,
		ylim : ylim,
		zlim : zlim,
		x_resolution:x_resolution,
		y_resolution:y_resolution

	};
};


//克里金生成矢量等值面
kriging.getVectorContour = function (gridInfo, breaks) {
	//像素坐标系的等值面
	var _contours = d3_contours()
		.size([gridInfo.n, gridInfo.m])
		.thresholds(breaks)
		(gridInfo.grid);
	//像素坐标系换算地理坐标系
	let dataset = {
		"type" : "FeatureCollection",
		"features" : []
	};
	var geoX_width=gridInfo.xlim[1]-gridInfo.xlim[0];
	var geoY_width=gridInfo.ylim[1]-gridInfo.ylim[0];
	_contours.forEach(contour => {
			contour.coordinates.forEach(polygon => {
						//polygon分内环和外环
						let _polygon = polygon.map(ring => {
									let _ring = ring.map(function (coor) {
										//像素坐标转地理坐标
										let lon = gridInfo.xlim[0] + geoX_width * (coor[0]*1.0 / gridInfo.n);
										let lat = gridInfo.ylim[0] + geoY_width * (coor[1]*1.0 / gridInfo.m);
										return [lon,lat];
									});
									return _ring;
								});
						dataset.features.push({
							"type" : "Feature",
							"properties" : {
								"contour_value" : contour.value
							},
							"geometry" : {
								"type" : "Polygon",
								"coordinates" : _polygon
							}
						});
			});
		});
	return dataset;
};
//克里金生成canvas图像
 kriging.drawCanvasContour = function(gridInfo,canvas,xlim,ylim,colors) {
	//清空画布
	var ctx = canvas.getContext("2d");
	ctx.clearRect(0, 0, canvas.width, canvas.height);
	//开始边界
	var range = [xlim[1]-xlim[0], ylim[1]-ylim[0], gridInfo.zlim[1]-gridInfo.zlim[0]];
	var n = gridInfo.n;
	var m = gridInfo.m;
	//根据分辨率，计算每个色块的宽高
	var wx = Math.ceil(gridInfo.x_resolution*canvas.width/(xlim[1]-xlim[0]));
	var wy = Math.ceil(gridInfo.y_resolution*canvas.height/(ylim[1]-ylim[0]));

	for(let i=0;i<m;i++)
	    for(let j=0;j<n;j++) {
			let _index=i*n+j;
			if(gridInfo.grid[_index]==undefined) 
				continue;
			let x = canvas.width*(j*gridInfo.x_resolution+gridInfo.xlim[0]-xlim[0])/range[0];
			let y = canvas.height*(1-(i*gridInfo.y_resolution+gridInfo.ylim[0]-ylim[0])/range[1]);
			let z = (gridInfo.grid[_index]-gridInfo.zlim[0])/range[2];
			if(z<0.0) 
				z = 0.0;
			else if(z>1.0) 
				z = 1.0;
			ctx.fillStyle = colors[Math.floor((colors.length-1)*z)];
			ctx.fillRect(Math.round(x-wx/2), Math.round(y-wy/2), wx, wy);
	    }
	};

var commonjsGlobal = typeof globalThis !== 'undefined' ? globalThis : typeof window !== 'undefined' ? window : typeof global !== 'undefined' ? global : typeof self !== 'undefined' ? self : {};

function unwrapExports (x) {
	return x && x.__esModule && Object.prototype.hasOwnProperty.call(x, 'default') ? x['default'] : x;
}

function createCommonjsModule(fn, module) {
	return module = { exports: {} }, fn(module, module.exports), module.exports;
}

var helpers = createCommonjsModule(function (module, exports) {
Object.defineProperty(exports, "__esModule", { value: true });
/**
 * @module helpers
 */
/**
 * Earth Radius used with the Harvesine formula and approximates using a spherical (non-ellipsoid) Earth.
 *
 * @memberof helpers
 * @type {number}
 */
exports.earthRadius = 6371008.8;
/**
 * Unit of measurement factors using a spherical (non-ellipsoid) earth radius.
 *
 * @memberof helpers
 * @type {Object}
 */
exports.factors = {
    centimeters: exports.earthRadius * 100,
    centimetres: exports.earthRadius * 100,
    degrees: exports.earthRadius / 111325,
    feet: exports.earthRadius * 3.28084,
    inches: exports.earthRadius * 39.370,
    kilometers: exports.earthRadius / 1000,
    kilometres: exports.earthRadius / 1000,
    meters: exports.earthRadius,
    metres: exports.earthRadius,
    miles: exports.earthRadius / 1609.344,
    millimeters: exports.earthRadius * 1000,
    millimetres: exports.earthRadius * 1000,
    nauticalmiles: exports.earthRadius / 1852,
    radians: 1,
    yards: exports.earthRadius / 1.0936,
};
/**
 * Units of measurement factors based on 1 meter.
 *
 * @memberof helpers
 * @type {Object}
 */
exports.unitsFactors = {
    centimeters: 100,
    centimetres: 100,
    degrees: 1 / 111325,
    feet: 3.28084,
    inches: 39.370,
    kilometers: 1 / 1000,
    kilometres: 1 / 1000,
    meters: 1,
    metres: 1,
    miles: 1 / 1609.344,
    millimeters: 1000,
    millimetres: 1000,
    nauticalmiles: 1 / 1852,
    radians: 1 / exports.earthRadius,
    yards: 1 / 1.0936,
};
/**
 * Area of measurement factors based on 1 square meter.
 *
 * @memberof helpers
 * @type {Object}
 */
exports.areaFactors = {
    acres: 0.000247105,
    centimeters: 10000,
    centimetres: 10000,
    feet: 10.763910417,
    inches: 1550.003100006,
    kilometers: 0.000001,
    kilometres: 0.000001,
    meters: 1,
    metres: 1,
    miles: 3.86e-7,
    millimeters: 1000000,
    millimetres: 1000000,
    yards: 1.195990046,
};
/**
 * Wraps a GeoJSON {@link Geometry} in a GeoJSON {@link Feature}.
 *
 * @name feature
 * @param {Geometry} geometry input geometry
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature} a GeoJSON Feature
 * @example
 * var geometry = {
 *   "type": "Point",
 *   "coordinates": [110, 50]
 * };
 *
 * var feature = turf.feature(geometry);
 *
 * //=feature
 */
function feature(geom, properties, options) {
    if (options === void 0) { options = {}; }
    var feat = { type: "Feature" };
    if (options.id === 0 || options.id) {
        feat.id = options.id;
    }
    if (options.bbox) {
        feat.bbox = options.bbox;
    }
    feat.properties = properties || {};
    feat.geometry = geom;
    return feat;
}
exports.feature = feature;
/**
 * Creates a GeoJSON {@link Geometry} from a Geometry string type & coordinates.
 * For GeometryCollection type use `helpers.geometryCollection`
 *
 * @name geometry
 * @param {string} type Geometry Type
 * @param {Array<any>} coordinates Coordinates
 * @param {Object} [options={}] Optional Parameters
 * @returns {Geometry} a GeoJSON Geometry
 * @example
 * var type = "Point";
 * var coordinates = [110, 50];
 * var geometry = turf.geometry(type, coordinates);
 * // => geometry
 */
function geometry(type, coordinates, options) {
    switch (type) {
        case "Point": return point(coordinates).geometry;
        case "LineString": return lineString(coordinates).geometry;
        case "Polygon": return polygon(coordinates).geometry;
        case "MultiPoint": return multiPoint(coordinates).geometry;
        case "MultiLineString": return multiLineString(coordinates).geometry;
        case "MultiPolygon": return multiPolygon(coordinates).geometry;
        default: throw new Error(type + " is invalid");
    }
}
exports.geometry = geometry;
/**
 * Creates a {@link Point} {@link Feature} from a Position.
 *
 * @name point
 * @param {Array<number>} coordinates longitude, latitude position (each in decimal degrees)
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<Point>} a Point feature
 * @example
 * var point = turf.point([-75.343, 39.984]);
 *
 * //=point
 */
function point(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    var geom = {
        type: "Point",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
exports.point = point;
/**
 * Creates a {@link Point} {@link FeatureCollection} from an Array of Point coordinates.
 *
 * @name points
 * @param {Array<Array<number>>} coordinates an array of Points
 * @param {Object} [properties={}] Translate these properties to each Feature
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north]
 * associated with the FeatureCollection
 * @param {string|number} [options.id] Identifier associated with the FeatureCollection
 * @returns {FeatureCollection<Point>} Point Feature
 * @example
 * var points = turf.points([
 *   [-75, 39],
 *   [-80, 45],
 *   [-78, 50]
 * ]);
 *
 * //=points
 */
function points(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    return featureCollection(coordinates.map(function (coords) {
        return point(coords, properties);
    }), options);
}
exports.points = points;
/**
 * Creates a {@link Polygon} {@link Feature} from an Array of LinearRings.
 *
 * @name polygon
 * @param {Array<Array<Array<number>>>} coordinates an array of LinearRings
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<Polygon>} Polygon Feature
 * @example
 * var polygon = turf.polygon([[[-5, 52], [-4, 56], [-2, 51], [-7, 54], [-5, 52]]], { name: 'poly1' });
 *
 * //=polygon
 */
function polygon(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    for (var _i = 0, coordinates_1 = coordinates; _i < coordinates_1.length; _i++) {
        var ring = coordinates_1[_i];
        if (ring.length < 4) {
            throw new Error("Each LinearRing of a Polygon must have 4 or more Positions.");
        }
        for (var j = 0; j < ring[ring.length - 1].length; j++) {
            // Check if first point of Polygon contains two numbers
            if (ring[ring.length - 1][j] !== ring[0][j]) {
                throw new Error("First and last Position are not equivalent.");
            }
        }
    }
    var geom = {
        type: "Polygon",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
exports.polygon = polygon;
/**
 * Creates a {@link Polygon} {@link FeatureCollection} from an Array of Polygon coordinates.
 *
 * @name polygons
 * @param {Array<Array<Array<Array<number>>>>} coordinates an array of Polygon coordinates
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the FeatureCollection
 * @returns {FeatureCollection<Polygon>} Polygon FeatureCollection
 * @example
 * var polygons = turf.polygons([
 *   [[[-5, 52], [-4, 56], [-2, 51], [-7, 54], [-5, 52]]],
 *   [[[-15, 42], [-14, 46], [-12, 41], [-17, 44], [-15, 42]]],
 * ]);
 *
 * //=polygons
 */
function polygons(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    return featureCollection(coordinates.map(function (coords) {
        return polygon(coords, properties);
    }), options);
}
exports.polygons = polygons;
/**
 * Creates a {@link LineString} {@link Feature} from an Array of Positions.
 *
 * @name lineString
 * @param {Array<Array<number>>} coordinates an array of Positions
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<LineString>} LineString Feature
 * @example
 * var linestring1 = turf.lineString([[-24, 63], [-23, 60], [-25, 65], [-20, 69]], {name: 'line 1'});
 * var linestring2 = turf.lineString([[-14, 43], [-13, 40], [-15, 45], [-10, 49]], {name: 'line 2'});
 *
 * //=linestring1
 * //=linestring2
 */
function lineString(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    if (coordinates.length < 2) {
        throw new Error("coordinates must be an array of two or more positions");
    }
    var geom = {
        type: "LineString",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
exports.lineString = lineString;
/**
 * Creates a {@link LineString} {@link FeatureCollection} from an Array of LineString coordinates.
 *
 * @name lineStrings
 * @param {Array<Array<Array<number>>>} coordinates an array of LinearRings
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north]
 * associated with the FeatureCollection
 * @param {string|number} [options.id] Identifier associated with the FeatureCollection
 * @returns {FeatureCollection<LineString>} LineString FeatureCollection
 * @example
 * var linestrings = turf.lineStrings([
 *   [[-24, 63], [-23, 60], [-25, 65], [-20, 69]],
 *   [[-14, 43], [-13, 40], [-15, 45], [-10, 49]]
 * ]);
 *
 * //=linestrings
 */
function lineStrings(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    return featureCollection(coordinates.map(function (coords) {
        return lineString(coords, properties);
    }), options);
}
exports.lineStrings = lineStrings;
/**
 * Takes one or more {@link Feature|Features} and creates a {@link FeatureCollection}.
 *
 * @name featureCollection
 * @param {Feature[]} features input features
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {FeatureCollection} FeatureCollection of Features
 * @example
 * var locationA = turf.point([-75.343, 39.984], {name: 'Location A'});
 * var locationB = turf.point([-75.833, 39.284], {name: 'Location B'});
 * var locationC = turf.point([-75.534, 39.123], {name: 'Location C'});
 *
 * var collection = turf.featureCollection([
 *   locationA,
 *   locationB,
 *   locationC
 * ]);
 *
 * //=collection
 */
function featureCollection(features, options) {
    if (options === void 0) { options = {}; }
    var fc = { type: "FeatureCollection" };
    if (options.id) {
        fc.id = options.id;
    }
    if (options.bbox) {
        fc.bbox = options.bbox;
    }
    fc.features = features;
    return fc;
}
exports.featureCollection = featureCollection;
/**
 * Creates a {@link Feature<MultiLineString>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name multiLineString
 * @param {Array<Array<Array<number>>>} coordinates an array of LineStrings
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<MultiLineString>} a MultiLineString feature
 * @throws {Error} if no coordinates are passed
 * @example
 * var multiLine = turf.multiLineString([[[0,0],[10,10]]]);
 *
 * //=multiLine
 */
function multiLineString(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    var geom = {
        type: "MultiLineString",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
exports.multiLineString = multiLineString;
/**
 * Creates a {@link Feature<MultiPoint>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name multiPoint
 * @param {Array<Array<number>>} coordinates an array of Positions
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<MultiPoint>} a MultiPoint feature
 * @throws {Error} if no coordinates are passed
 * @example
 * var multiPt = turf.multiPoint([[0,0],[10,10]]);
 *
 * //=multiPt
 */
function multiPoint(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    var geom = {
        type: "MultiPoint",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
exports.multiPoint = multiPoint;
/**
 * Creates a {@link Feature<MultiPolygon>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name multiPolygon
 * @param {Array<Array<Array<Array<number>>>>} coordinates an array of Polygons
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<MultiPolygon>} a multipolygon feature
 * @throws {Error} if no coordinates are passed
 * @example
 * var multiPoly = turf.multiPolygon([[[[0,0],[0,10],[10,10],[10,0],[0,0]]]]);
 *
 * //=multiPoly
 *
 */
function multiPolygon(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    var geom = {
        type: "MultiPolygon",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
exports.multiPolygon = multiPolygon;
/**
 * Creates a {@link Feature<GeometryCollection>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name geometryCollection
 * @param {Array<Geometry>} geometries an array of GeoJSON Geometries
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<GeometryCollection>} a GeoJSON GeometryCollection Feature
 * @example
 * var pt = turf.geometry("Point", [100, 0]);
 * var line = turf.geometry("LineString", [[101, 0], [102, 1]]);
 * var collection = turf.geometryCollection([pt, line]);
 *
 * // => collection
 */
function geometryCollection(geometries, properties, options) {
    if (options === void 0) { options = {}; }
    var geom = {
        type: "GeometryCollection",
        geometries: geometries,
    };
    return feature(geom, properties, options);
}
exports.geometryCollection = geometryCollection;
/**
 * Round number to precision
 *
 * @param {number} num Number
 * @param {number} [precision=0] Precision
 * @returns {number} rounded number
 * @example
 * turf.round(120.4321)
 * //=120
 *
 * turf.round(120.4321, 2)
 * //=120.43
 */
function round(num, precision) {
    if (precision === void 0) { precision = 0; }
    if (precision && !(precision >= 0)) {
        throw new Error("precision must be a positive number");
    }
    var multiplier = Math.pow(10, precision || 0);
    return Math.round(num * multiplier) / multiplier;
}
exports.round = round;
/**
 * Convert a distance measurement (assuming a spherical Earth) from radians to a more friendly unit.
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
 *
 * @name radiansToLength
 * @param {number} radians in radians across the sphere
 * @param {string} [units="kilometers"] can be degrees, radians, miles, or kilometers inches, yards, metres,
 * meters, kilometres, kilometers.
 * @returns {number} distance
 */
function radiansToLength(radians, units) {
    if (units === void 0) { units = "kilometers"; }
    var factor = exports.factors[units];
    if (!factor) {
        throw new Error(units + " units is invalid");
    }
    return radians * factor;
}
exports.radiansToLength = radiansToLength;
/**
 * Convert a distance measurement (assuming a spherical Earth) from a real-world unit into radians
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
 *
 * @name lengthToRadians
 * @param {number} distance in real units
 * @param {string} [units="kilometers"] can be degrees, radians, miles, or kilometers inches, yards, metres,
 * meters, kilometres, kilometers.
 * @returns {number} radians
 */
function lengthToRadians(distance, units) {
    if (units === void 0) { units = "kilometers"; }
    var factor = exports.factors[units];
    if (!factor) {
        throw new Error(units + " units is invalid");
    }
    return distance / factor;
}
exports.lengthToRadians = lengthToRadians;
/**
 * Convert a distance measurement (assuming a spherical Earth) from a real-world unit into degrees
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, centimeters, kilometres, feet
 *
 * @name lengthToDegrees
 * @param {number} distance in real units
 * @param {string} [units="kilometers"] can be degrees, radians, miles, or kilometers inches, yards, metres,
 * meters, kilometres, kilometers.
 * @returns {number} degrees
 */
function lengthToDegrees(distance, units) {
    return radiansToDegrees(lengthToRadians(distance, units));
}
exports.lengthToDegrees = lengthToDegrees;
/**
 * Converts any bearing angle from the north line direction (positive clockwise)
 * and returns an angle between 0-360 degrees (positive clockwise), 0 being the north line
 *
 * @name bearingToAzimuth
 * @param {number} bearing angle, between -180 and +180 degrees
 * @returns {number} angle between 0 and 360 degrees
 */
function bearingToAzimuth(bearing) {
    var angle = bearing % 360;
    if (angle < 0) {
        angle += 360;
    }
    return angle;
}
exports.bearingToAzimuth = bearingToAzimuth;
/**
 * Converts an angle in radians to degrees
 *
 * @name radiansToDegrees
 * @param {number} radians angle in radians
 * @returns {number} degrees between 0 and 360 degrees
 */
function radiansToDegrees(radians) {
    var degrees = radians % (2 * Math.PI);
    return degrees * 180 / Math.PI;
}
exports.radiansToDegrees = radiansToDegrees;
/**
 * Converts an angle in degrees to radians
 *
 * @name degreesToRadians
 * @param {number} degrees angle between 0 and 360 degrees
 * @returns {number} angle in radians
 */
function degreesToRadians(degrees) {
    var radians = degrees % 360;
    return radians * Math.PI / 180;
}
exports.degreesToRadians = degreesToRadians;
/**
 * Converts a length to the requested unit.
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
 *
 * @param {number} length to be converted
 * @param {Units} [originalUnit="kilometers"] of the length
 * @param {Units} [finalUnit="kilometers"] returned unit
 * @returns {number} the converted length
 */
function convertLength(length, originalUnit, finalUnit) {
    if (originalUnit === void 0) { originalUnit = "kilometers"; }
    if (finalUnit === void 0) { finalUnit = "kilometers"; }
    if (!(length >= 0)) {
        throw new Error("length must be a positive number");
    }
    return radiansToLength(lengthToRadians(length, originalUnit), finalUnit);
}
exports.convertLength = convertLength;
/**
 * Converts a area to the requested unit.
 * Valid units: kilometers, kilometres, meters, metres, centimetres, millimeters, acres, miles, yards, feet, inches
 * @param {number} area to be converted
 * @param {Units} [originalUnit="meters"] of the distance
 * @param {Units} [finalUnit="kilometers"] returned unit
 * @returns {number} the converted distance
 */
function convertArea(area, originalUnit, finalUnit) {
    if (originalUnit === void 0) { originalUnit = "meters"; }
    if (finalUnit === void 0) { finalUnit = "kilometers"; }
    if (!(area >= 0)) {
        throw new Error("area must be a positive number");
    }
    var startFactor = exports.areaFactors[originalUnit];
    if (!startFactor) {
        throw new Error("invalid original units");
    }
    var finalFactor = exports.areaFactors[finalUnit];
    if (!finalFactor) {
        throw new Error("invalid final units");
    }
    return (area / startFactor) * finalFactor;
}
exports.convertArea = convertArea;
/**
 * isNumber
 *
 * @param {*} num Number to validate
 * @returns {boolean} true/false
 * @example
 * turf.isNumber(123)
 * //=true
 * turf.isNumber('foo')
 * //=false
 */
function isNumber(num) {
    return !isNaN(num) && num !== null && !Array.isArray(num) && !/^\s*$/.test(num);
}
exports.isNumber = isNumber;
/**
 * isObject
 *
 * @param {*} input variable to validate
 * @returns {boolean} true/false
 * @example
 * turf.isObject({elevation: 10})
 * //=true
 * turf.isObject('foo')
 * //=false
 */
function isObject(input) {
    return (!!input) && (input.constructor === Object);
}
exports.isObject = isObject;
/**
 * Validate BBox
 *
 * @private
 * @param {Array<number>} bbox BBox to validate
 * @returns {void}
 * @throws Error if BBox is not valid
 * @example
 * validateBBox([-180, -40, 110, 50])
 * //=OK
 * validateBBox([-180, -40])
 * //=Error
 * validateBBox('Foo')
 * //=Error
 * validateBBox(5)
 * //=Error
 * validateBBox(null)
 * //=Error
 * validateBBox(undefined)
 * //=Error
 */
function validateBBox(bbox) {
    if (!bbox) {
        throw new Error("bbox is required");
    }
    if (!Array.isArray(bbox)) {
        throw new Error("bbox must be an Array");
    }
    if (bbox.length !== 4 && bbox.length !== 6) {
        throw new Error("bbox must be an Array of 4 or 6 numbers");
    }
    bbox.forEach(function (num) {
        if (!isNumber(num)) {
            throw new Error("bbox must only contain numbers");
        }
    });
}
exports.validateBBox = validateBBox;
/**
 * Validate Id
 *
 * @private
 * @param {string|number} id Id to validate
 * @returns {void}
 * @throws Error if Id is not valid
 * @example
 * validateId([-180, -40, 110, 50])
 * //=Error
 * validateId([-180, -40])
 * //=Error
 * validateId('Foo')
 * //=OK
 * validateId(5)
 * //=OK
 * validateId(null)
 * //=Error
 * validateId(undefined)
 * //=Error
 */
function validateId(id) {
    if (!id) {
        throw new Error("id is required");
    }
    if (["string", "number"].indexOf(typeof id) === -1) {
        throw new Error("id must be a number or a string");
    }
}
exports.validateId = validateId;
// Deprecated methods
function radians2degrees() {
    throw new Error("method has been renamed to `radiansToDegrees`");
}
exports.radians2degrees = radians2degrees;
function degrees2radians() {
    throw new Error("method has been renamed to `degreesToRadians`");
}
exports.degrees2radians = degrees2radians;
function distanceToDegrees() {
    throw new Error("method has been renamed to `lengthToDegrees`");
}
exports.distanceToDegrees = distanceToDegrees;
function distanceToRadians() {
    throw new Error("method has been renamed to `lengthToRadians`");
}
exports.distanceToRadians = distanceToRadians;
function radiansToDistance() {
    throw new Error("method has been renamed to `radiansToLength`");
}
exports.radiansToDistance = radiansToDistance;
function bearingToAngle() {
    throw new Error("method has been renamed to `bearingToAzimuth`");
}
exports.bearingToAngle = bearingToAngle;
function convertDistance() {
    throw new Error("method has been renamed to `convertLength`");
}
exports.convertDistance = convertDistance;
});

unwrapExports(helpers);
var helpers_1 = helpers.earthRadius;
var helpers_2 = helpers.factors;
var helpers_3 = helpers.unitsFactors;
var helpers_4 = helpers.areaFactors;
var helpers_5 = helpers.feature;
var helpers_6 = helpers.geometry;
var helpers_7 = helpers.point;
var helpers_8 = helpers.points;
var helpers_9 = helpers.polygon;
var helpers_10 = helpers.polygons;
var helpers_11 = helpers.lineString;
var helpers_12 = helpers.lineStrings;
var helpers_13 = helpers.featureCollection;
var helpers_14 = helpers.multiLineString;
var helpers_15 = helpers.multiPoint;
var helpers_16 = helpers.multiPolygon;
var helpers_17 = helpers.geometryCollection;
var helpers_18 = helpers.round;
var helpers_19 = helpers.radiansToLength;
var helpers_20 = helpers.lengthToRadians;
var helpers_21 = helpers.lengthToDegrees;
var helpers_22 = helpers.bearingToAzimuth;
var helpers_23 = helpers.radiansToDegrees;
var helpers_24 = helpers.degreesToRadians;
var helpers_25 = helpers.convertLength;
var helpers_26 = helpers.convertArea;
var helpers_27 = helpers.isNumber;
var helpers_28 = helpers.isObject;
var helpers_29 = helpers.validateBBox;
var helpers_30 = helpers.validateId;
var helpers_31 = helpers.radians2degrees;
var helpers_32 = helpers.degrees2radians;
var helpers_33 = helpers.distanceToDegrees;
var helpers_34 = helpers.distanceToRadians;
var helpers_35 = helpers.radiansToDistance;
var helpers_36 = helpers.bearingToAngle;
var helpers_37 = helpers.convertDistance;

var helpers$1 = createCommonjsModule(function (module, exports) {
Object.defineProperty(exports, "__esModule", { value: true });
/**
 * @module helpers
 */
/**
 * Earth Radius used with the Harvesine formula and approximates using a spherical (non-ellipsoid) Earth.
 *
 * @memberof helpers
 * @type {number}
 */
exports.earthRadius = 6371008.8;
/**
 * Unit of measurement factors using a spherical (non-ellipsoid) earth radius.
 *
 * @memberof helpers
 * @type {Object}
 */
exports.factors = {
    centimeters: exports.earthRadius * 100,
    centimetres: exports.earthRadius * 100,
    degrees: exports.earthRadius / 111325,
    feet: exports.earthRadius * 3.28084,
    inches: exports.earthRadius * 39.370,
    kilometers: exports.earthRadius / 1000,
    kilometres: exports.earthRadius / 1000,
    meters: exports.earthRadius,
    metres: exports.earthRadius,
    miles: exports.earthRadius / 1609.344,
    millimeters: exports.earthRadius * 1000,
    millimetres: exports.earthRadius * 1000,
    nauticalmiles: exports.earthRadius / 1852,
    radians: 1,
    yards: exports.earthRadius / 1.0936,
};
/**
 * Units of measurement factors based on 1 meter.
 *
 * @memberof helpers
 * @type {Object}
 */
exports.unitsFactors = {
    centimeters: 100,
    centimetres: 100,
    degrees: 1 / 111325,
    feet: 3.28084,
    inches: 39.370,
    kilometers: 1 / 1000,
    kilometres: 1 / 1000,
    meters: 1,
    metres: 1,
    miles: 1 / 1609.344,
    millimeters: 1000,
    millimetres: 1000,
    nauticalmiles: 1 / 1852,
    radians: 1 / exports.earthRadius,
    yards: 1 / 1.0936,
};
/**
 * Area of measurement factors based on 1 square meter.
 *
 * @memberof helpers
 * @type {Object}
 */
exports.areaFactors = {
    acres: 0.000247105,
    centimeters: 10000,
    centimetres: 10000,
    feet: 10.763910417,
    inches: 1550.003100006,
    kilometers: 0.000001,
    kilometres: 0.000001,
    meters: 1,
    metres: 1,
    miles: 3.86e-7,
    millimeters: 1000000,
    millimetres: 1000000,
    yards: 1.195990046,
};
/**
 * Wraps a GeoJSON {@link Geometry} in a GeoJSON {@link Feature}.
 *
 * @name feature
 * @param {Geometry} geometry input geometry
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature} a GeoJSON Feature
 * @example
 * var geometry = {
 *   "type": "Point",
 *   "coordinates": [110, 50]
 * };
 *
 * var feature = turf.feature(geometry);
 *
 * //=feature
 */
function feature(geom, properties, options) {
    if (options === void 0) { options = {}; }
    var feat = { type: "Feature" };
    if (options.id === 0 || options.id) {
        feat.id = options.id;
    }
    if (options.bbox) {
        feat.bbox = options.bbox;
    }
    feat.properties = properties || {};
    feat.geometry = geom;
    return feat;
}
exports.feature = feature;
/**
 * Creates a GeoJSON {@link Geometry} from a Geometry string type & coordinates.
 * For GeometryCollection type use `helpers.geometryCollection`
 *
 * @name geometry
 * @param {string} type Geometry Type
 * @param {Array<any>} coordinates Coordinates
 * @param {Object} [options={}] Optional Parameters
 * @returns {Geometry} a GeoJSON Geometry
 * @example
 * var type = "Point";
 * var coordinates = [110, 50];
 * var geometry = turf.geometry(type, coordinates);
 * // => geometry
 */
function geometry(type, coordinates, options) {
    switch (type) {
        case "Point": return point(coordinates).geometry;
        case "LineString": return lineString(coordinates).geometry;
        case "Polygon": return polygon(coordinates).geometry;
        case "MultiPoint": return multiPoint(coordinates).geometry;
        case "MultiLineString": return multiLineString(coordinates).geometry;
        case "MultiPolygon": return multiPolygon(coordinates).geometry;
        default: throw new Error(type + " is invalid");
    }
}
exports.geometry = geometry;
/**
 * Creates a {@link Point} {@link Feature} from a Position.
 *
 * @name point
 * @param {Array<number>} coordinates longitude, latitude position (each in decimal degrees)
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<Point>} a Point feature
 * @example
 * var point = turf.point([-75.343, 39.984]);
 *
 * //=point
 */
function point(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    var geom = {
        type: "Point",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
exports.point = point;
/**
 * Creates a {@link Point} {@link FeatureCollection} from an Array of Point coordinates.
 *
 * @name points
 * @param {Array<Array<number>>} coordinates an array of Points
 * @param {Object} [properties={}] Translate these properties to each Feature
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north]
 * associated with the FeatureCollection
 * @param {string|number} [options.id] Identifier associated with the FeatureCollection
 * @returns {FeatureCollection<Point>} Point Feature
 * @example
 * var points = turf.points([
 *   [-75, 39],
 *   [-80, 45],
 *   [-78, 50]
 * ]);
 *
 * //=points
 */
function points(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    return featureCollection(coordinates.map(function (coords) {
        return point(coords, properties);
    }), options);
}
exports.points = points;
/**
 * Creates a {@link Polygon} {@link Feature} from an Array of LinearRings.
 *
 * @name polygon
 * @param {Array<Array<Array<number>>>} coordinates an array of LinearRings
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<Polygon>} Polygon Feature
 * @example
 * var polygon = turf.polygon([[[-5, 52], [-4, 56], [-2, 51], [-7, 54], [-5, 52]]], { name: 'poly1' });
 *
 * //=polygon
 */
function polygon(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    for (var _i = 0, coordinates_1 = coordinates; _i < coordinates_1.length; _i++) {
        var ring = coordinates_1[_i];
        if (ring.length < 4) {
            throw new Error("Each LinearRing of a Polygon must have 4 or more Positions.");
        }
        for (var j = 0; j < ring[ring.length - 1].length; j++) {
            // Check if first point of Polygon contains two numbers
            if (ring[ring.length - 1][j] !== ring[0][j]) {
                throw new Error("First and last Position are not equivalent.");
            }
        }
    }
    var geom = {
        type: "Polygon",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
exports.polygon = polygon;
/**
 * Creates a {@link Polygon} {@link FeatureCollection} from an Array of Polygon coordinates.
 *
 * @name polygons
 * @param {Array<Array<Array<Array<number>>>>} coordinates an array of Polygon coordinates
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the FeatureCollection
 * @returns {FeatureCollection<Polygon>} Polygon FeatureCollection
 * @example
 * var polygons = turf.polygons([
 *   [[[-5, 52], [-4, 56], [-2, 51], [-7, 54], [-5, 52]]],
 *   [[[-15, 42], [-14, 46], [-12, 41], [-17, 44], [-15, 42]]],
 * ]);
 *
 * //=polygons
 */
function polygons(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    return featureCollection(coordinates.map(function (coords) {
        return polygon(coords, properties);
    }), options);
}
exports.polygons = polygons;
/**
 * Creates a {@link LineString} {@link Feature} from an Array of Positions.
 *
 * @name lineString
 * @param {Array<Array<number>>} coordinates an array of Positions
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<LineString>} LineString Feature
 * @example
 * var linestring1 = turf.lineString([[-24, 63], [-23, 60], [-25, 65], [-20, 69]], {name: 'line 1'});
 * var linestring2 = turf.lineString([[-14, 43], [-13, 40], [-15, 45], [-10, 49]], {name: 'line 2'});
 *
 * //=linestring1
 * //=linestring2
 */
function lineString(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    if (coordinates.length < 2) {
        throw new Error("coordinates must be an array of two or more positions");
    }
    var geom = {
        type: "LineString",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
exports.lineString = lineString;
/**
 * Creates a {@link LineString} {@link FeatureCollection} from an Array of LineString coordinates.
 *
 * @name lineStrings
 * @param {Array<Array<Array<number>>>} coordinates an array of LinearRings
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north]
 * associated with the FeatureCollection
 * @param {string|number} [options.id] Identifier associated with the FeatureCollection
 * @returns {FeatureCollection<LineString>} LineString FeatureCollection
 * @example
 * var linestrings = turf.lineStrings([
 *   [[-24, 63], [-23, 60], [-25, 65], [-20, 69]],
 *   [[-14, 43], [-13, 40], [-15, 45], [-10, 49]]
 * ]);
 *
 * //=linestrings
 */
function lineStrings(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    return featureCollection(coordinates.map(function (coords) {
        return lineString(coords, properties);
    }), options);
}
exports.lineStrings = lineStrings;
/**
 * Takes one or more {@link Feature|Features} and creates a {@link FeatureCollection}.
 *
 * @name featureCollection
 * @param {Feature[]} features input features
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {FeatureCollection} FeatureCollection of Features
 * @example
 * var locationA = turf.point([-75.343, 39.984], {name: 'Location A'});
 * var locationB = turf.point([-75.833, 39.284], {name: 'Location B'});
 * var locationC = turf.point([-75.534, 39.123], {name: 'Location C'});
 *
 * var collection = turf.featureCollection([
 *   locationA,
 *   locationB,
 *   locationC
 * ]);
 *
 * //=collection
 */
function featureCollection(features, options) {
    if (options === void 0) { options = {}; }
    var fc = { type: "FeatureCollection" };
    if (options.id) {
        fc.id = options.id;
    }
    if (options.bbox) {
        fc.bbox = options.bbox;
    }
    fc.features = features;
    return fc;
}
exports.featureCollection = featureCollection;
/**
 * Creates a {@link Feature<MultiLineString>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name multiLineString
 * @param {Array<Array<Array<number>>>} coordinates an array of LineStrings
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<MultiLineString>} a MultiLineString feature
 * @throws {Error} if no coordinates are passed
 * @example
 * var multiLine = turf.multiLineString([[[0,0],[10,10]]]);
 *
 * //=multiLine
 */
function multiLineString(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    var geom = {
        type: "MultiLineString",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
exports.multiLineString = multiLineString;
/**
 * Creates a {@link Feature<MultiPoint>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name multiPoint
 * @param {Array<Array<number>>} coordinates an array of Positions
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<MultiPoint>} a MultiPoint feature
 * @throws {Error} if no coordinates are passed
 * @example
 * var multiPt = turf.multiPoint([[0,0],[10,10]]);
 *
 * //=multiPt
 */
function multiPoint(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    var geom = {
        type: "MultiPoint",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
exports.multiPoint = multiPoint;
/**
 * Creates a {@link Feature<MultiPolygon>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name multiPolygon
 * @param {Array<Array<Array<Array<number>>>>} coordinates an array of Polygons
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<MultiPolygon>} a multipolygon feature
 * @throws {Error} if no coordinates are passed
 * @example
 * var multiPoly = turf.multiPolygon([[[[0,0],[0,10],[10,10],[10,0],[0,0]]]]);
 *
 * //=multiPoly
 *
 */
function multiPolygon(coordinates, properties, options) {
    if (options === void 0) { options = {}; }
    var geom = {
        type: "MultiPolygon",
        coordinates: coordinates,
    };
    return feature(geom, properties, options);
}
exports.multiPolygon = multiPolygon;
/**
 * Creates a {@link Feature<GeometryCollection>} based on a
 * coordinate array. Properties can be added optionally.
 *
 * @name geometryCollection
 * @param {Array<Geometry>} geometries an array of GeoJSON Geometries
 * @param {Object} [properties={}] an Object of key-value pairs to add as properties
 * @param {Object} [options={}] Optional Parameters
 * @param {Array<number>} [options.bbox] Bounding Box Array [west, south, east, north] associated with the Feature
 * @param {string|number} [options.id] Identifier associated with the Feature
 * @returns {Feature<GeometryCollection>} a GeoJSON GeometryCollection Feature
 * @example
 * var pt = turf.geometry("Point", [100, 0]);
 * var line = turf.geometry("LineString", [[101, 0], [102, 1]]);
 * var collection = turf.geometryCollection([pt, line]);
 *
 * // => collection
 */
function geometryCollection(geometries, properties, options) {
    if (options === void 0) { options = {}; }
    var geom = {
        type: "GeometryCollection",
        geometries: geometries,
    };
    return feature(geom, properties, options);
}
exports.geometryCollection = geometryCollection;
/**
 * Round number to precision
 *
 * @param {number} num Number
 * @param {number} [precision=0] Precision
 * @returns {number} rounded number
 * @example
 * turf.round(120.4321)
 * //=120
 *
 * turf.round(120.4321, 2)
 * //=120.43
 */
function round(num, precision) {
    if (precision === void 0) { precision = 0; }
    if (precision && !(precision >= 0)) {
        throw new Error("precision must be a positive number");
    }
    var multiplier = Math.pow(10, precision || 0);
    return Math.round(num * multiplier) / multiplier;
}
exports.round = round;
/**
 * Convert a distance measurement (assuming a spherical Earth) from radians to a more friendly unit.
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
 *
 * @name radiansToLength
 * @param {number} radians in radians across the sphere
 * @param {string} [units="kilometers"] can be degrees, radians, miles, or kilometers inches, yards, metres,
 * meters, kilometres, kilometers.
 * @returns {number} distance
 */
function radiansToLength(radians, units) {
    if (units === void 0) { units = "kilometers"; }
    var factor = exports.factors[units];
    if (!factor) {
        throw new Error(units + " units is invalid");
    }
    return radians * factor;
}
exports.radiansToLength = radiansToLength;
/**
 * Convert a distance measurement (assuming a spherical Earth) from a real-world unit into radians
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
 *
 * @name lengthToRadians
 * @param {number} distance in real units
 * @param {string} [units="kilometers"] can be degrees, radians, miles, or kilometers inches, yards, metres,
 * meters, kilometres, kilometers.
 * @returns {number} radians
 */
function lengthToRadians(distance, units) {
    if (units === void 0) { units = "kilometers"; }
    var factor = exports.factors[units];
    if (!factor) {
        throw new Error(units + " units is invalid");
    }
    return distance / factor;
}
exports.lengthToRadians = lengthToRadians;
/**
 * Convert a distance measurement (assuming a spherical Earth) from a real-world unit into degrees
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, centimeters, kilometres, feet
 *
 * @name lengthToDegrees
 * @param {number} distance in real units
 * @param {string} [units="kilometers"] can be degrees, radians, miles, or kilometers inches, yards, metres,
 * meters, kilometres, kilometers.
 * @returns {number} degrees
 */
function lengthToDegrees(distance, units) {
    return radiansToDegrees(lengthToRadians(distance, units));
}
exports.lengthToDegrees = lengthToDegrees;
/**
 * Converts any bearing angle from the north line direction (positive clockwise)
 * and returns an angle between 0-360 degrees (positive clockwise), 0 being the north line
 *
 * @name bearingToAzimuth
 * @param {number} bearing angle, between -180 and +180 degrees
 * @returns {number} angle between 0 and 360 degrees
 */
function bearingToAzimuth(bearing) {
    var angle = bearing % 360;
    if (angle < 0) {
        angle += 360;
    }
    return angle;
}
exports.bearingToAzimuth = bearingToAzimuth;
/**
 * Converts an angle in radians to degrees
 *
 * @name radiansToDegrees
 * @param {number} radians angle in radians
 * @returns {number} degrees between 0 and 360 degrees
 */
function radiansToDegrees(radians) {
    var degrees = radians % (2 * Math.PI);
    return degrees * 180 / Math.PI;
}
exports.radiansToDegrees = radiansToDegrees;
/**
 * Converts an angle in degrees to radians
 *
 * @name degreesToRadians
 * @param {number} degrees angle between 0 and 360 degrees
 * @returns {number} angle in radians
 */
function degreesToRadians(degrees) {
    var radians = degrees % 360;
    return radians * Math.PI / 180;
}
exports.degreesToRadians = degreesToRadians;
/**
 * Converts a length to the requested unit.
 * Valid units: miles, nauticalmiles, inches, yards, meters, metres, kilometers, centimeters, feet
 *
 * @param {number} length to be converted
 * @param {Units} [originalUnit="kilometers"] of the length
 * @param {Units} [finalUnit="kilometers"] returned unit
 * @returns {number} the converted length
 */
function convertLength(length, originalUnit, finalUnit) {
    if (originalUnit === void 0) { originalUnit = "kilometers"; }
    if (finalUnit === void 0) { finalUnit = "kilometers"; }
    if (!(length >= 0)) {
        throw new Error("length must be a positive number");
    }
    return radiansToLength(lengthToRadians(length, originalUnit), finalUnit);
}
exports.convertLength = convertLength;
/**
 * Converts a area to the requested unit.
 * Valid units: kilometers, kilometres, meters, metres, centimetres, millimeters, acres, miles, yards, feet, inches
 * @param {number} area to be converted
 * @param {Units} [originalUnit="meters"] of the distance
 * @param {Units} [finalUnit="kilometers"] returned unit
 * @returns {number} the converted distance
 */
function convertArea(area, originalUnit, finalUnit) {
    if (originalUnit === void 0) { originalUnit = "meters"; }
    if (finalUnit === void 0) { finalUnit = "kilometers"; }
    if (!(area >= 0)) {
        throw new Error("area must be a positive number");
    }
    var startFactor = exports.areaFactors[originalUnit];
    if (!startFactor) {
        throw new Error("invalid original units");
    }
    var finalFactor = exports.areaFactors[finalUnit];
    if (!finalFactor) {
        throw new Error("invalid final units");
    }
    return (area / startFactor) * finalFactor;
}
exports.convertArea = convertArea;
/**
 * isNumber
 *
 * @param {*} num Number to validate
 * @returns {boolean} true/false
 * @example
 * turf.isNumber(123)
 * //=true
 * turf.isNumber('foo')
 * //=false
 */
function isNumber(num) {
    return !isNaN(num) && num !== null && !Array.isArray(num) && !/^\s*$/.test(num);
}
exports.isNumber = isNumber;
/**
 * isObject
 *
 * @param {*} input variable to validate
 * @returns {boolean} true/false
 * @example
 * turf.isObject({elevation: 10})
 * //=true
 * turf.isObject('foo')
 * //=false
 */
function isObject(input) {
    return (!!input) && (input.constructor === Object);
}
exports.isObject = isObject;
/**
 * Validate BBox
 *
 * @private
 * @param {Array<number>} bbox BBox to validate
 * @returns {void}
 * @throws Error if BBox is not valid
 * @example
 * validateBBox([-180, -40, 110, 50])
 * //=OK
 * validateBBox([-180, -40])
 * //=Error
 * validateBBox('Foo')
 * //=Error
 * validateBBox(5)
 * //=Error
 * validateBBox(null)
 * //=Error
 * validateBBox(undefined)
 * //=Error
 */
function validateBBox(bbox) {
    if (!bbox) {
        throw new Error("bbox is required");
    }
    if (!Array.isArray(bbox)) {
        throw new Error("bbox must be an Array");
    }
    if (bbox.length !== 4 && bbox.length !== 6) {
        throw new Error("bbox must be an Array of 4 or 6 numbers");
    }
    bbox.forEach(function (num) {
        if (!isNumber(num)) {
            throw new Error("bbox must only contain numbers");
        }
    });
}
exports.validateBBox = validateBBox;
/**
 * Validate Id
 *
 * @private
 * @param {string|number} id Id to validate
 * @returns {void}
 * @throws Error if Id is not valid
 * @example
 * validateId([-180, -40, 110, 50])
 * //=Error
 * validateId([-180, -40])
 * //=Error
 * validateId('Foo')
 * //=OK
 * validateId(5)
 * //=OK
 * validateId(null)
 * //=Error
 * validateId(undefined)
 * //=Error
 */
function validateId(id) {
    if (!id) {
        throw new Error("id is required");
    }
    if (["string", "number"].indexOf(typeof id) === -1) {
        throw new Error("id must be a number or a string");
    }
}
exports.validateId = validateId;
// Deprecated methods
function radians2degrees() {
    throw new Error("method has been renamed to `radiansToDegrees`");
}
exports.radians2degrees = radians2degrees;
function degrees2radians() {
    throw new Error("method has been renamed to `degreesToRadians`");
}
exports.degrees2radians = degrees2radians;
function distanceToDegrees() {
    throw new Error("method has been renamed to `lengthToDegrees`");
}
exports.distanceToDegrees = distanceToDegrees;
function distanceToRadians() {
    throw new Error("method has been renamed to `lengthToRadians`");
}
exports.distanceToRadians = distanceToRadians;
function radiansToDistance() {
    throw new Error("method has been renamed to `radiansToLength`");
}
exports.radiansToDistance = radiansToDistance;
function bearingToAngle() {
    throw new Error("method has been renamed to `bearingToAzimuth`");
}
exports.bearingToAngle = bearingToAngle;
function convertDistance() {
    throw new Error("method has been renamed to `convertLength`");
}
exports.convertDistance = convertDistance;
});

unwrapExports(helpers$1);
var helpers_1$1 = helpers$1.earthRadius;
var helpers_2$1 = helpers$1.factors;
var helpers_3$1 = helpers$1.unitsFactors;
var helpers_4$1 = helpers$1.areaFactors;
var helpers_5$1 = helpers$1.feature;
var helpers_6$1 = helpers$1.geometry;
var helpers_7$1 = helpers$1.point;
var helpers_8$1 = helpers$1.points;
var helpers_9$1 = helpers$1.polygon;
var helpers_10$1 = helpers$1.polygons;
var helpers_11$1 = helpers$1.lineString;
var helpers_12$1 = helpers$1.lineStrings;
var helpers_13$1 = helpers$1.featureCollection;
var helpers_14$1 = helpers$1.multiLineString;
var helpers_15$1 = helpers$1.multiPoint;
var helpers_16$1 = helpers$1.multiPolygon;
var helpers_17$1 = helpers$1.geometryCollection;
var helpers_18$1 = helpers$1.round;
var helpers_19$1 = helpers$1.radiansToLength;
var helpers_20$1 = helpers$1.lengthToRadians;
var helpers_21$1 = helpers$1.lengthToDegrees;
var helpers_22$1 = helpers$1.bearingToAzimuth;
var helpers_23$1 = helpers$1.radiansToDegrees;
var helpers_24$1 = helpers$1.degreesToRadians;
var helpers_25$1 = helpers$1.convertLength;
var helpers_26$1 = helpers$1.convertArea;
var helpers_27$1 = helpers$1.isNumber;
var helpers_28$1 = helpers$1.isObject;
var helpers_29$1 = helpers$1.validateBBox;
var helpers_30$1 = helpers$1.validateId;
var helpers_31$1 = helpers$1.radians2degrees;
var helpers_32$1 = helpers$1.degrees2radians;
var helpers_33$1 = helpers$1.distanceToDegrees;
var helpers_34$1 = helpers$1.distanceToRadians;
var helpers_35$1 = helpers$1.radiansToDistance;
var helpers_36$1 = helpers$1.bearingToAngle;
var helpers_37$1 = helpers$1.convertDistance;

var invariant = createCommonjsModule(function (module, exports) {
Object.defineProperty(exports, "__esModule", { value: true });

/**
 * Unwrap a coordinate from a Point Feature, Geometry or a single coordinate.
 *
 * @name getCoord
 * @param {Array<number>|Geometry<Point>|Feature<Point>} coord GeoJSON Point or an Array of numbers
 * @returns {Array<number>} coordinates
 * @example
 * var pt = turf.point([10, 10]);
 *
 * var coord = turf.getCoord(pt);
 * //= [10, 10]
 */
function getCoord(coord) {
    if (!coord) {
        throw new Error("coord is required");
    }
    if (!Array.isArray(coord)) {
        if (coord.type === "Feature" && coord.geometry !== null && coord.geometry.type === "Point") {
            return coord.geometry.coordinates;
        }
        if (coord.type === "Point") {
            return coord.coordinates;
        }
    }
    if (Array.isArray(coord) && coord.length >= 2 && !Array.isArray(coord[0]) && !Array.isArray(coord[1])) {
        return coord;
    }
    throw new Error("coord must be GeoJSON Point or an Array of numbers");
}
exports.getCoord = getCoord;
/**
 * Unwrap coordinates from a Feature, Geometry Object or an Array
 *
 * @name getCoords
 * @param {Array<any>|Geometry|Feature} coords Feature, Geometry Object or an Array
 * @returns {Array<any>} coordinates
 * @example
 * var poly = turf.polygon([[[119.32, -8.7], [119.55, -8.69], [119.51, -8.54], [119.32, -8.7]]]);
 *
 * var coords = turf.getCoords(poly);
 * //= [[[119.32, -8.7], [119.55, -8.69], [119.51, -8.54], [119.32, -8.7]]]
 */
function getCoords(coords) {
    if (Array.isArray(coords)) {
        return coords;
    }
    // Feature
    if (coords.type === "Feature") {
        if (coords.geometry !== null) {
            return coords.geometry.coordinates;
        }
    }
    else {
        // Geometry
        if (coords.coordinates) {
            return coords.coordinates;
        }
    }
    throw new Error("coords must be GeoJSON Feature, Geometry Object or an Array");
}
exports.getCoords = getCoords;
/**
 * Checks if coordinates contains a number
 *
 * @name containsNumber
 * @param {Array<any>} coordinates GeoJSON Coordinates
 * @returns {boolean} true if Array contains a number
 */
function containsNumber(coordinates) {
    if (coordinates.length > 1 && helpers$1.isNumber(coordinates[0]) && helpers$1.isNumber(coordinates[1])) {
        return true;
    }
    if (Array.isArray(coordinates[0]) && coordinates[0].length) {
        return containsNumber(coordinates[0]);
    }
    throw new Error("coordinates must only contain numbers");
}
exports.containsNumber = containsNumber;
/**
 * Enforce expectations about types of GeoJSON objects for Turf.
 *
 * @name geojsonType
 * @param {GeoJSON} value any GeoJSON object
 * @param {string} type expected GeoJSON type
 * @param {string} name name of calling function
 * @throws {Error} if value is not the expected type.
 */
function geojsonType(value, type, name) {
    if (!type || !name) {
        throw new Error("type and name required");
    }
    if (!value || value.type !== type) {
        throw new Error("Invalid input to " + name + ": must be a " + type + ", given " + value.type);
    }
}
exports.geojsonType = geojsonType;
/**
 * Enforce expectations about types of {@link Feature} inputs for Turf.
 * Internally this uses {@link geojsonType} to judge geometry types.
 *
 * @name featureOf
 * @param {Feature} feature a feature with an expected geometry type
 * @param {string} type expected GeoJSON type
 * @param {string} name name of calling function
 * @throws {Error} error if value is not the expected type.
 */
function featureOf(feature, type, name) {
    if (!feature) {
        throw new Error("No feature passed");
    }
    if (!name) {
        throw new Error(".featureOf() requires a name");
    }
    if (!feature || feature.type !== "Feature" || !feature.geometry) {
        throw new Error("Invalid input to " + name + ", Feature with geometry required");
    }
    if (!feature.geometry || feature.geometry.type !== type) {
        throw new Error("Invalid input to " + name + ": must be a " + type + ", given " + feature.geometry.type);
    }
}
exports.featureOf = featureOf;
/**
 * Enforce expectations about types of {@link FeatureCollection} inputs for Turf.
 * Internally this uses {@link geojsonType} to judge geometry types.
 *
 * @name collectionOf
 * @param {FeatureCollection} featureCollection a FeatureCollection for which features will be judged
 * @param {string} type expected GeoJSON type
 * @param {string} name name of calling function
 * @throws {Error} if value is not the expected type.
 */
function collectionOf(featureCollection, type, name) {
    if (!featureCollection) {
        throw new Error("No featureCollection passed");
    }
    if (!name) {
        throw new Error(".collectionOf() requires a name");
    }
    if (!featureCollection || featureCollection.type !== "FeatureCollection") {
        throw new Error("Invalid input to " + name + ", FeatureCollection required");
    }
    for (var _i = 0, _a = featureCollection.features; _i < _a.length; _i++) {
        var feature = _a[_i];
        if (!feature || feature.type !== "Feature" || !feature.geometry) {
            throw new Error("Invalid input to " + name + ", Feature with geometry required");
        }
        if (!feature.geometry || feature.geometry.type !== type) {
            throw new Error("Invalid input to " + name + ": must be a " + type + ", given " + feature.geometry.type);
        }
    }
}
exports.collectionOf = collectionOf;
/**
 * Get Geometry from Feature or Geometry Object
 *
 * @param {Feature|Geometry} geojson GeoJSON Feature or Geometry Object
 * @returns {Geometry|null} GeoJSON Geometry Object
 * @throws {Error} if geojson is not a Feature or Geometry Object
 * @example
 * var point = {
 *   "type": "Feature",
 *   "properties": {},
 *   "geometry": {
 *     "type": "Point",
 *     "coordinates": [110, 40]
 *   }
 * }
 * var geom = turf.getGeom(point)
 * //={"type": "Point", "coordinates": [110, 40]}
 */
function getGeom(geojson) {
    if (geojson.type === "Feature") {
        return geojson.geometry;
    }
    return geojson;
}
exports.getGeom = getGeom;
/**
 * Get GeoJSON object's type, Geometry type is prioritize.
 *
 * @param {GeoJSON} geojson GeoJSON object
 * @param {string} [name="geojson"] name of the variable to display in error message
 * @returns {string} GeoJSON type
 * @example
 * var point = {
 *   "type": "Feature",
 *   "properties": {},
 *   "geometry": {
 *     "type": "Point",
 *     "coordinates": [110, 40]
 *   }
 * }
 * var geom = turf.getType(point)
 * //="Point"
 */
function getType(geojson, name) {
    if (geojson.type === "FeatureCollection") {
        return "FeatureCollection";
    }
    if (geojson.type === "GeometryCollection") {
        return "GeometryCollection";
    }
    if (geojson.type === "Feature" && geojson.geometry !== null) {
        return geojson.geometry.type;
    }
    return geojson.type;
}
exports.getType = getType;
});

unwrapExports(invariant);
var invariant_1 = invariant.getCoord;
var invariant_2 = invariant.getCoords;
var invariant_3 = invariant.containsNumber;
var invariant_4 = invariant.geojsonType;
var invariant_5 = invariant.featureOf;
var invariant_6 = invariant.collectionOf;
var invariant_7 = invariant.getGeom;
var invariant_8 = invariant.getType;

function DEFAULT_COMPARE (a, b) { return a > b ? 1 : a < b ? -1 : 0; }

class SplayTree {

  constructor(compare = DEFAULT_COMPARE, noDuplicates = false) {
    this._compare = compare;
    this._root = null;
    this._size = 0;
    this._noDuplicates = !!noDuplicates;
  }


  rotateLeft(x) {
    var y = x.right;
    if (y) {
      x.right = y.left;
      if (y.left) y.left.parent = x;
      y.parent = x.parent;
    }

    if (!x.parent)                this._root = y;
    else if (x === x.parent.left) x.parent.left = y;
    else                          x.parent.right = y;
    if (y) y.left = x;
    x.parent = y;
  }


  rotateRight(x) {
    var y = x.left;
    if (y) {
      x.left = y.right;
      if (y.right) y.right.parent = x;
      y.parent = x.parent;
    }

    if (!x.parent)               this._root = y;
    else if(x === x.parent.left) x.parent.left = y;
    else                         x.parent.right = y;
    if (y) y.right = x;
    x.parent = y;
  }


  _splay(x) {
    while (x.parent) {
      var p = x.parent;
      if (!p.parent) {
        if (p.left === x) this.rotateRight(p);
        else              this.rotateLeft(p);
      } else if (p.left === x && p.parent.left === p) {
        this.rotateRight(p.parent);
        this.rotateRight(p);
      } else if (p.right === x && p.parent.right === p) {
        this.rotateLeft(p.parent);
        this.rotateLeft(p);
      } else if (p.left === x && p.parent.right === p) {
        this.rotateRight(p);
        this.rotateLeft(p);
      } else {
        this.rotateLeft(p);
        this.rotateRight(p);
      }
    }
  }


  splay(x) {
    var p, gp, ggp, l, r;

    while (x.parent) {
      p = x.parent;
      gp = p.parent;

      if (gp && gp.parent) {
        ggp = gp.parent;
        if (ggp.left === gp) ggp.left  = x;
        else                 ggp.right = x;
        x.parent = ggp;
      } else {
        x.parent = null;
        this._root = x;
      }

      l = x.left; r = x.right;

      if (x === p.left) { // left
        if (gp) {
          if (gp.left === p) {
            /* zig-zig */
            if (p.right) {
              gp.left = p.right;
              gp.left.parent = gp;
            } else gp.left = null;

            p.right   = gp;
            gp.parent = p;
          } else {
            /* zig-zag */
            if (l) {
              gp.right = l;
              l.parent = gp;
            } else gp.right = null;

            x.left    = gp;
            gp.parent = x;
          }
        }
        if (r) {
          p.left = r;
          r.parent = p;
        } else p.left = null;

        x.right  = p;
        p.parent = x;
      } else { // right
        if (gp) {
          if (gp.right === p) {
            /* zig-zig */
            if (p.left) {
              gp.right = p.left;
              gp.right.parent = gp;
            } else gp.right = null;

            p.left = gp;
            gp.parent = p;
          } else {
            /* zig-zag */
            if (r) {
              gp.left = r;
              r.parent = gp;
            } else gp.left = null;

            x.right   = gp;
            gp.parent = x;
          }
        }
        if (l) {
          p.right = l;
          l.parent = p;
        } else p.right = null;

        x.left   = p;
        p.parent = x;
      }
    }
  }


  replace(u, v) {
    if (!u.parent) this._root = v;
    else if (u === u.parent.left) u.parent.left = v;
    else u.parent.right = v;
    if (v) v.parent = u.parent;
  }


  minNode(u = this._root) {
    if (u) while (u.left) u = u.left;
    return u;
  }


  maxNode(u = this._root) {
    if (u) while (u.right) u = u.right;
    return u;
  }


  insert(key, data) {
    var z = this._root;
    var p = null;
    var comp = this._compare;
    var cmp;

    if (this._noDuplicates) {
      while (z) {
        p = z;
        cmp = comp(z.key, key);
        if (cmp === 0) return;
        else if (comp(z.key, key) < 0) z = z.right;
        else z = z.left;
      }
    } else {
      while (z) {
        p = z;
        if (comp(z.key, key) < 0) z = z.right;
        else z = z.left;
      }
    }

    z = { key, data, left: null, right: null, parent: p };

    if (!p)                          this._root = z;
    else if (comp(p.key, z.key) < 0) p.right = z;
    else                             p.left  = z;

    this.splay(z);
    this._size++;
    return z;
  }


  find (key) {
    var z    = this._root;
    var comp = this._compare;
    while (z) {
      var cmp = comp(z.key, key);
      if      (cmp < 0) z = z.right;
      else if (cmp > 0) z = z.left;
      else              return z;
    }
    return null;
  }

  /**
   * Whether the tree contains a node with the given key
   * @param  {Key} key
   * @return {boolean} true/false
   */
  contains (key) {
    var node       = this._root;
    var comparator = this._compare;
    while (node)  {
      var cmp = comparator(key, node.key);
      if      (cmp === 0) return true;
      else if (cmp < 0)   node = node.left;
      else                node = node.right;
    }

    return false;
  }


  remove (key) {
    var z = this.find(key);

    if (!z) return false;

    this.splay(z);

    if (!z.left) this.replace(z, z.right);
    else if (!z.right) this.replace(z, z.left);
    else {
      var y = this.minNode(z.right);
      if (y.parent !== z) {
        this.replace(y, y.right);
        y.right = z.right;
        y.right.parent = y;
      }
      this.replace(z, y);
      y.left = z.left;
      y.left.parent = y;
    }

    this._size--;
    return true;
  }


  removeNode(z) {
    if (!z) return false;

    this.splay(z);

    if (!z.left) this.replace(z, z.right);
    else if (!z.right) this.replace(z, z.left);
    else {
      var y = this.minNode(z.right);
      if (y.parent !== z) {
        this.replace(y, y.right);
        y.right = z.right;
        y.right.parent = y;
      }
      this.replace(z, y);
      y.left = z.left;
      y.left.parent = y;
    }

    this._size--;
    return true;
  }


  erase (key) {
    var z = this.find(key);
    if (!z) return;

    this.splay(z);

    var s = z.left;
    var t = z.right;

    var sMax = null;
    if (s) {
      s.parent = null;
      sMax = this.maxNode(s);
      this.splay(sMax);
      this._root = sMax;
    }
    if (t) {
      if (s) sMax.right = t;
      else   this._root = t;
      t.parent = sMax;
    }

    this._size--;
  }

  /**
   * Removes and returns the node with smallest key
   * @return {?Node}
   */
  pop () {
    var node = this._root, returnValue = null;
    if (node) {
      while (node.left) node = node.left;
      returnValue = { key: node.key, data: node.data };
      this.remove(node.key);
    }
    return returnValue;
  }


  /* eslint-disable class-methods-use-this */

  /**
   * Successor node
   * @param  {Node} node
   * @return {?Node}
   */
  next (node) {
    var successor = node;
    if (successor) {
      if (successor.right) {
        successor = successor.right;
        while (successor && successor.left) successor = successor.left;
      } else {
        successor = node.parent;
        while (successor && successor.right === node) {
          node = successor; successor = successor.parent;
        }
      }
    }
    return successor;
  }


  /**
   * Predecessor node
   * @param  {Node} node
   * @return {?Node}
   */
  prev (node) {
    var predecessor = node;
    if (predecessor) {
      if (predecessor.left) {
        predecessor = predecessor.left;
        while (predecessor && predecessor.right) predecessor = predecessor.right;
      } else {
        predecessor = node.parent;
        while (predecessor && predecessor.left === node) {
          node = predecessor;
          predecessor = predecessor.parent;
        }
      }
    }
    return predecessor;
  }
  /* eslint-enable class-methods-use-this */


  /**
   * @param  {forEachCallback} callback
   * @return {SplayTree}
   */
  forEach(callback) {
    var current = this._root;
    var s = [], done = false, i = 0;

    while (!done) {
      // Reach the left most Node of the current Node
      if (current) {
        // Place pointer to a tree node on the stack
        // before traversing the node's left subtree
        s.push(current);
        current = current.left;
      } else {
        // BackTrack from the empty subtree and visit the Node
        // at the top of the stack; however, if the stack is
        // empty you are done
        if (s.length > 0) {
          current = s.pop();
          callback(current, i++);

          // We have visited the node and its left
          // subtree. Now, it's right subtree's turn
          current = current.right;
        } else done = true;
      }
    }
    return this;
  }


  /**
   * Walk key range from `low` to `high`. Stops if `fn` returns a value.
   * @param  {Key}      low
   * @param  {Key}      high
   * @param  {Function} fn
   * @param  {*?}       ctx
   * @return {SplayTree}
   */
  range(low, high, fn, ctx) {
    const Q = [];
    const compare = this._compare;
    let node = this._root, cmp;

    while (Q.length !== 0 || node) {
      if (node) {
        Q.push(node);
        node = node.left;
      } else {
        node = Q.pop();
        cmp = compare(node.key, high);
        if (cmp > 0) {
          break;
        } else if (compare(node.key, low) >= 0) {
          if (fn.call(ctx, node)) return this; // stop if smth is returned
        }
        node = node.right;
      }
    }
    return this;
  }

  /**
   * Returns all keys in order
   * @return {Array<Key>}
   */
  keys () {
    var current = this._root;
    var s = [], r = [], done = false;

    while (!done) {
      if (current) {
        s.push(current);
        current = current.left;
      } else {
        if (s.length > 0) {
          current = s.pop();
          r.push(current.key);
          current = current.right;
        } else done = true;
      }
    }
    return r;
  }


  /**
   * Returns `data` fields of all nodes in order.
   * @return {Array<Value>}
   */
  values () {
    var current = this._root;
    var s = [], r = [], done = false;

    while (!done) {
      if (current) {
        s.push(current);
        current = current.left;
      } else {
        if (s.length > 0) {
          current = s.pop();
          r.push(current.data);
          current = current.right;
        } else done = true;
      }
    }
    return r;
  }


  /**
   * Returns node at given index
   * @param  {number} index
   * @return {?Node}
   */
  at (index) {
    // removed after a consideration, more misleading than useful
    // index = index % this.size;
    // if (index < 0) index = this.size - index;

    var current = this._root;
    var s = [], done = false, i = 0;

    while (!done) {
      if (current) {
        s.push(current);
        current = current.left;
      } else {
        if (s.length > 0) {
          current = s.pop();
          if (i === index) return current;
          i++;
          current = current.right;
        } else done = true;
      }
    }
    return null;
  }

  /**
   * Bulk-load items. Both array have to be same size
   * @param  {Array<Key>}    keys
   * @param  {Array<Value>}  [values]
   * @param  {Boolean}       [presort=false] Pre-sort keys and values, using
   *                                         tree's comparator. Sorting is done
   *                                         in-place
   * @return {AVLTree}
   */
  load(keys = [], values = [], presort = false) {
    if (this._size !== 0) throw new Error('bulk-load: tree is not empty');
    const size = keys.length;
    if (presort) sort(keys, values, 0, size - 1, this._compare);
    this._root = loadRecursive(null, keys, values, 0, size);
    this._size = size;
    return this;
  }


  min() {
    var node = this.minNode(this._root);
    if (node) return node.key;
    else      return null;
  }


  max() {
    var node = this.maxNode(this._root);
    if (node) return node.key;
    else      return null;
  }

  isEmpty() { return this._root === null; }
  get size() { return this._size; }


  /**
   * Create a tree and load it with items
   * @param  {Array<Key>}          keys
   * @param  {Array<Value>?}        [values]

   * @param  {Function?}            [comparator]
   * @param  {Boolean?}             [presort=false] Pre-sort keys and values, using
   *                                               tree's comparator. Sorting is done
   *                                               in-place
   * @param  {Boolean?}             [noDuplicates=false]   Allow duplicates
   * @return {SplayTree}
   */
  static createTree(keys, values, comparator, presort, noDuplicates) {
    return new SplayTree(comparator, noDuplicates).load(keys, values, presort);
  }
}


function loadRecursive (parent, keys, values, start, end) {
  const size = end - start;
  if (size > 0) {
    const middle = start + Math.floor(size / 2);
    const key    = keys[middle];
    const data   = values[middle];
    const node   = { key, data, parent };
    node.left    = loadRecursive(node, keys, values, start, middle);
    node.right   = loadRecursive(node, keys, values, middle + 1, end);
    return node;
  }
  return null;
}


function sort(keys, values, left, right, compare) {
  if (left >= right) return;

  const pivot = keys[(left + right) >> 1];
  let i = left - 1;
  let j = right + 1;

  while (true) {
    do i++; while (compare(keys[i], pivot) < 0);
    do j--; while (compare(keys[j], pivot) > 0);
    if (i >= j) break;

    let tmp = keys[i];
    keys[i] = keys[j];
    keys[j] = tmp;

    tmp = values[i];
    values[i] = values[j];
    values[j] = tmp;
  }

  sort(keys, values,  left,     j, compare);
  sort(keys, values, j + 1, right, compare);
}

const NORMAL               = 0;
const NON_CONTRIBUTING     = 1;
const SAME_TRANSITION      = 2;
const DIFFERENT_TRANSITION = 3;

const INTERSECTION = 0;
const UNION        = 1;
const DIFFERENCE   = 2;
const XOR          = 3;

/**
 * @param  {SweepEvent} event
 * @param  {SweepEvent} prev
 * @param  {Operation} operation
 */
function computeFields (event, prev, operation) {
  // compute inOut and otherInOut fields
  if (prev === null) {
    event.inOut      = false;
    event.otherInOut = true;

  // previous line segment in sweepline belongs to the same polygon
  } else {
    if (event.isSubject === prev.isSubject) {
      event.inOut      = !prev.inOut;
      event.otherInOut = prev.otherInOut;

    // previous line segment in sweepline belongs to the clipping polygon
    } else {
      event.inOut      = !prev.otherInOut;
      event.otherInOut = prev.isVertical() ? !prev.inOut : prev.inOut;
    }

    // compute prevInResult field
    if (prev) {
      event.prevInResult = (!inResult(prev, operation) || prev.isVertical())
        ? prev.prevInResult : prev;
    }
  }

  // check if the line segment belongs to the Boolean operation
  event.inResult = inResult(event, operation);
}


/* eslint-disable indent */
function inResult(event, operation) {
  switch (event.type) {
    case NORMAL:
      switch (operation) {
        case INTERSECTION:
          return !event.otherInOut;
        case UNION:
          return event.otherInOut;
        case DIFFERENCE:
          // return (event.isSubject && !event.otherInOut) ||
          //         (!event.isSubject && event.otherInOut);
          return (event.isSubject && event.otherInOut) ||
                  (!event.isSubject && !event.otherInOut);
        case XOR:
          return true;
      }
      break;
    case SAME_TRANSITION:
      return operation === INTERSECTION || operation === UNION;
    case DIFFERENT_TRANSITION:
      return operation === DIFFERENCE;
    case NON_CONTRIBUTING:
      return false;
  }
  return false;
}
/* eslint-enable indent */

class SweepEvent {


  /**
   * Sweepline event
   *
   * @class {SweepEvent}
   * @param {Array.<Number>}  point
   * @param {Boolean}         left
   * @param {SweepEvent=}     otherEvent
   * @param {Boolean}         isSubject
   * @param {Number}          edgeType
   */
  constructor (point, left, otherEvent, isSubject, edgeType) {

    /**
     * Is left endpoint?
     * @type {Boolean}
     */
    this.left = left;

    /**
     * @type {Array.<Number>}
     */
    this.point = point;

    /**
     * Other edge reference
     * @type {SweepEvent}
     */
    this.otherEvent = otherEvent;

    /**
     * Belongs to source or clipping polygon
     * @type {Boolean}
     */
    this.isSubject = isSubject;

    /**
     * Edge contribution type
     * @type {Number}
     */
    this.type = edgeType || NORMAL;


    /**
     * In-out transition for the sweepline crossing polygon
     * @type {Boolean}
     */
    this.inOut = false;


    /**
     * @type {Boolean}
     */
    this.otherInOut = false;

    /**
     * Previous event in result?
     * @type {SweepEvent}
     */
    this.prevInResult = null;

    /**
     * Does event belong to result?
     * @type {Boolean}
     */
    this.inResult = false;


    // connection step

    /**
     * @type {Boolean}
     */
    this.resultInOut = false;

    this.isExteriorRing = true;
  }


  /**
   * @param  {Array.<Number>}  p
   * @return {Boolean}
   */
  isBelow (p) {
    const p0 = this.point, p1 = this.otherEvent.point;
    return this.left
      ? (p0[0] - p[0]) * (p1[1] - p[1]) - (p1[0] - p[0]) * (p0[1] - p[1]) > 0
      // signedArea(this.point, this.otherEvent.point, p) > 0 :
      : (p1[0] - p[0]) * (p0[1] - p[1]) - (p0[0] - p[0]) * (p1[1] - p[1]) > 0;
      //signedArea(this.otherEvent.point, this.point, p) > 0;
  }


  /**
   * @param  {Array.<Number>}  p
   * @return {Boolean}
   */
  isAbove (p) {
    return !this.isBelow(p);
  }


  /**
   * @return {Boolean}
   */
  isVertical () {
    return this.point[0] === this.otherEvent.point[0];
  }


  clone () {
    const copy = new SweepEvent(
      this.point, this.left, this.otherEvent, this.isSubject, this.type);

    copy.inResult       = this.inResult;
    copy.prevInResult   = this.prevInResult;
    copy.isExteriorRing = this.isExteriorRing;
    copy.inOut          = this.inOut;
    copy.otherInOut     = this.otherInOut;

    return copy;
  }
}

function equals(p1, p2) {
  if (p1[0] === p2[0]) {
    if (p1[1] === p2[1]) {
      return true;
    } else {
      return false;
    }
  }
  return false;
}

// const EPSILON = 1e-9;
// const abs = Math.abs;
// TODO https://github.com/w8r/martinez/issues/6#issuecomment-262847164
// Precision problem.
//
// module.exports = function equals(p1, p2) {
//   return abs(p1[0] - p2[0]) <= EPSILON && abs(p1[1] - p2[1]) <= EPSILON;
// };

/**
 * Signed area of the triangle (p0, p1, p2)
 * @param  {Array.<Number>} p0
 * @param  {Array.<Number>} p1
 * @param  {Array.<Number>} p2
 * @return {Number}
 */
function signedArea(p0, p1, p2) {
  return (p0[0] - p2[0]) * (p1[1] - p2[1]) - (p1[0] - p2[0]) * (p0[1] - p2[1]);
}

/**
 * @param  {SweepEvent} e1
 * @param  {SweepEvent} e2
 * @return {Number}
 */
function compareEvents(e1, e2) {
  const p1 = e1.point;
  const p2 = e2.point;

  // Different x-coordinate
  if (p1[0] > p2[0]) return 1;
  if (p1[0] < p2[0]) return -1;

  // Different points, but same x-coordinate
  // Event with lower y-coordinate is processed first
  if (p1[1] !== p2[1]) return p1[1] > p2[1] ? 1 : -1;

  return specialCases(e1, e2, p1);
}


/* eslint-disable no-unused-vars */
function specialCases(e1, e2, p1, p2) {
  // Same coordinates, but one is a left endpoint and the other is
  // a right endpoint. The right endpoint is processed first
  if (e1.left !== e2.left)
    return e1.left ? 1 : -1;

  // const p2 = e1.otherEvent.point, p3 = e2.otherEvent.point;
  // const sa = (p1[0] - p3[0]) * (p2[1] - p3[1]) - (p2[0] - p3[0]) * (p1[1] - p3[1])
  // Same coordinates, both events
  // are left endpoints or right endpoints.
  // not collinear
  if (signedArea(p1, e1.otherEvent.point, e2.otherEvent.point) !== 0) {
    // the event associate to the bottom segment is processed first
    return (!e1.isBelow(e2.otherEvent.point)) ? 1 : -1;
  }

  return (!e1.isSubject && e2.isSubject) ? 1 : -1;
}
/* eslint-enable no-unused-vars */

/**
 * @param  {SweepEvent} se
 * @param  {Array.<Number>} p
 * @param  {Queue} queue
 * @return {Queue}
 */
function divideSegment(se, p, queue)  {
  const r = new SweepEvent(p, false, se,            se.isSubject);
  const l = new SweepEvent(p, true,  se.otherEvent, se.isSubject);

  /* eslint-disable no-console */
  if (equals(se.point, se.otherEvent.point)) {

    console.warn('what is that, a collapsed segment?', se);
  }
  /* eslint-enable no-console */

  r.contourId = l.contourId = se.contourId;

  // avoid a rounding error. The left event would be processed after the right event
  if (compareEvents(l, se.otherEvent) > 0) {
    se.otherEvent.left = true;
    l.left = false;
  }

  // avoid a rounding error. The left event would be processed after the right event
  // if (compareEvents(se, r) > 0) {}

  se.otherEvent.otherEvent = l;
  se.otherEvent = r;

  queue.push(l);
  queue.push(r);

  return queue;
}

//const EPS = 1e-9;

/**
 * Finds the magnitude of the cross product of two vectors (if we pretend
 * they're in three dimensions)
 *
 * @param {Object} a First vector
 * @param {Object} b Second vector
 * @private
 * @returns {Number} The magnitude of the cross product
 */
function crossProduct(a, b) {
  return (a[0] * b[1]) - (a[1] * b[0]);
}

/**
 * Finds the dot product of two vectors.
 *
 * @param {Object} a First vector
 * @param {Object} b Second vector
 * @private
 * @returns {Number} The dot product
 */
function dotProduct(a, b) {
  return (a[0] * b[0]) + (a[1] * b[1]);
}

/**
 * Finds the intersection (if any) between two line segments a and b, given the
 * line segments' end points a1, a2 and b1, b2.
 *
 * This algorithm is based on Schneider and Eberly.
 * http://www.cimec.org.ar/~ncalvo/Schneider_Eberly.pdf
 * Page 244.
 *
 * @param {Array.<Number>} a1 point of first line
 * @param {Array.<Number>} a2 point of first line
 * @param {Array.<Number>} b1 point of second line
 * @param {Array.<Number>} b2 point of second line
 * @param {Boolean=}       noEndpointTouch whether to skip single touchpoints
 *                                         (meaning connected segments) as
 *                                         intersections
 * @returns {Array.<Array.<Number>>|Null} If the lines intersect, the point of
 * intersection. If they overlap, the two end points of the overlapping segment.
 * Otherwise, null.
 */
function intersection (a1, a2, b1, b2, noEndpointTouch) {
  // The algorithm expects our lines in the form P + sd, where P is a point,
  // s is on the interval [0, 1], and d is a vector.
  // We are passed two points. P can be the first point of each pair. The
  // vector, then, could be thought of as the distance (in x and y components)
  // from the first point to the second point.
  // So first, let's make our vectors:
  const va = [a2[0] - a1[0], a2[1] - a1[1]];
  const vb = [b2[0] - b1[0], b2[1] - b1[1]];
  // We also define a function to convert back to regular point form:

  /* eslint-disable arrow-body-style */

  function toPoint(p, s, d) {
    return [
      p[0] + s * d[0],
      p[1] + s * d[1]
    ];
  }

  /* eslint-enable arrow-body-style */

  // The rest is pretty much a straight port of the algorithm.
  const e = [b1[0] - a1[0], b1[1] - a1[1]];
  let kross    = crossProduct(va, vb);
  let sqrKross = kross * kross;
  const sqrLenA  = dotProduct(va, va);
  //const sqrLenB  = dotProduct(vb, vb);

  // Check for line intersection. This works because of the properties of the
  // cross product -- specifically, two vectors are parallel if and only if the
  // cross product is the 0 vector. The full calculation involves relative error
  // to account for possible very small line segments. See Schneider & Eberly
  // for details.
  if (sqrKross > 0/* EPS * sqrLenB * sqLenA */) {
    // If they're not parallel, then (because these are line segments) they
    // still might not actually intersect. This code checks that the
    // intersection point of the lines is actually on both line segments.
    const s = crossProduct(e, vb) / kross;
    if (s < 0 || s > 1) {
      // not on line segment a
      return null;
    }
    const t = crossProduct(e, va) / kross;
    if (t < 0 || t > 1) {
      // not on line segment b
      return null;
    }
    if (s === 0 || s === 1) {
      // on an endpoint of line segment a
      return noEndpointTouch ? null : [toPoint(a1, s, va)];
    }
    if (t === 0 || t === 1) {
      // on an endpoint of line segment b
      return noEndpointTouch ? null : [toPoint(b1, t, vb)];
    }
    return [toPoint(a1, s, va)];
  }

  // If we've reached this point, then the lines are either parallel or the
  // same, but the segments could overlap partially or fully, or not at all.
  // So we need to find the overlap, if any. To do that, we can use e, which is
  // the (vector) difference between the two initial points. If this is parallel
  // with the line itself, then the two lines are the same line, and there will
  // be overlap.
  //const sqrLenE = dotProduct(e, e);
  kross = crossProduct(e, va);
  sqrKross = kross * kross;

  if (sqrKross > 0 /* EPS * sqLenB * sqLenE */) {
  // Lines are just parallel, not the same. No overlap.
    return null;
  }

  const sa = dotProduct(va, e) / sqrLenA;
  const sb = sa + dotProduct(va, vb) / sqrLenA;
  const smin = Math.min(sa, sb);
  const smax = Math.max(sa, sb);

  // this is, essentially, the FindIntersection acting on floats from
  // Schneider & Eberly, just inlined into this function.
  if (smin <= 1 && smax >= 0) {

    // overlap on an end point
    if (smin === 1) {
      return noEndpointTouch ? null : [toPoint(a1, smin > 0 ? smin : 0, va)];
    }

    if (smax === 0) {
      return noEndpointTouch ? null : [toPoint(a1, smax < 1 ? smax : 1, va)];
    }

    if (noEndpointTouch && smin === 0 && smax === 1) return null;

    // There's overlap on a segment -- two points of intersection. Return both.
    return [
      toPoint(a1, smin > 0 ? smin : 0, va),
      toPoint(a1, smax < 1 ? smax : 1, va)
    ];
  }

  return null;
}

/**
 * @param  {SweepEvent} se1
 * @param  {SweepEvent} se2
 * @param  {Queue}      queue
 * @return {Number}
 */
function possibleIntersection (se1, se2, queue) {
  // that disallows self-intersecting polygons,
  // did cost us half a day, so I'll leave it
  // out of respect
  // if (se1.isSubject === se2.isSubject) return;
  const inter = intersection(
    se1.point, se1.otherEvent.point,
    se2.point, se2.otherEvent.point
  );

  const nintersections = inter ? inter.length : 0;
  if (nintersections === 0) return 0; // no intersection

  // the line segments intersect at an endpoint of both line segments
  if ((nintersections === 1) &&
      (equals(se1.point, se2.point) ||
       equals(se1.otherEvent.point, se2.otherEvent.point))) {
    return 0;
  }

  if (nintersections === 2 && se1.isSubject === se2.isSubject) {
    // if(se1.contourId === se2.contourId){
    // console.warn('Edges of the same polygon overlap',
    //   se1.point, se1.otherEvent.point, se2.point, se2.otherEvent.point);
    // }
    //throw new Error('Edges of the same polygon overlap');
    return 0;
  }

  // The line segments associated to se1 and se2 intersect
  if (nintersections === 1) {

    // if the intersection point is not an endpoint of se1
    if (!equals(se1.point, inter[0]) && !equals(se1.otherEvent.point, inter[0])) {
      divideSegment(se1, inter[0], queue);
    }

    // if the intersection point is not an endpoint of se2
    if (!equals(se2.point, inter[0]) && !equals(se2.otherEvent.point, inter[0])) {
      divideSegment(se2, inter[0], queue);
    }
    return 1;
  }

  // The line segments associated to se1 and se2 overlap
  const events        = [];
  let leftCoincide  = false;
  let rightCoincide = false;

  if (equals(se1.point, se2.point)) {
    leftCoincide = true; // linked
  } else if (compareEvents(se1, se2) === 1) {
    events.push(se2, se1);
  } else {
    events.push(se1, se2);
  }

  if (equals(se1.otherEvent.point, se2.otherEvent.point)) {
    rightCoincide = true;
  } else if (compareEvents(se1.otherEvent, se2.otherEvent) === 1) {
    events.push(se2.otherEvent, se1.otherEvent);
  } else {
    events.push(se1.otherEvent, se2.otherEvent);
  }

  if ((leftCoincide && rightCoincide) || leftCoincide) {
    // both line segments are equal or share the left endpoint
    se2.type = NON_CONTRIBUTING;
    se1.type = (se2.inOut === se1.inOut)
      ? SAME_TRANSITION : DIFFERENT_TRANSITION;

    if (leftCoincide && !rightCoincide) {
      // honestly no idea, but changing events selection from [2, 1]
      // to [0, 1] fixes the overlapping self-intersecting polygons issue
      divideSegment(events[1].otherEvent, events[0].point, queue);
    }
    return 2;
  }

  // the line segments share the right endpoint
  if (rightCoincide) {
    divideSegment(events[0], events[1].point, queue);
    return 3;
  }

  // no line segment includes totally the other one
  if (events[0] !== events[3].otherEvent) {
    divideSegment(events[0], events[1].point, queue);
    divideSegment(events[1], events[2].point, queue);
    return 3;
  }

  // one line segment includes the other one
  divideSegment(events[0], events[1].point, queue);
  divideSegment(events[3].otherEvent, events[2].point, queue);

  return 3;
}

/**
 * @param  {SweepEvent} le1
 * @param  {SweepEvent} le2
 * @return {Number}
 */
function compareSegments(le1, le2) {
  if (le1 === le2) return 0;

  // Segments are not collinear
  if (signedArea(le1.point, le1.otherEvent.point, le2.point) !== 0 ||
    signedArea(le1.point, le1.otherEvent.point, le2.otherEvent.point) !== 0) {

    // If they share their left endpoint use the right endpoint to sort
    if (equals(le1.point, le2.point)) return le1.isBelow(le2.otherEvent.point) ? -1 : 1;

    // Different left endpoint: use the left endpoint to sort
    if (le1.point[0] === le2.point[0]) return le1.point[1] < le2.point[1] ? -1 : 1;

    // has the line segment associated to e1 been inserted
    // into S after the line segment associated to e2 ?
    if (compareEvents(le1, le2) === 1) return le2.isAbove(le1.point) ? -1 : 1;

    // The line segment associated to e2 has been inserted
    // into S after the line segment associated to e1
    return le1.isBelow(le2.point) ? -1 : 1;
  }

  if (le1.isSubject === le2.isSubject) { // same polygon
    let p1 = le1.point, p2 = le2.point;
    if (p1[0] === p2[0] && p1[1] === p2[1]/*equals(le1.point, le2.point)*/) {
      p1 = le1.otherEvent.point; p2 = le2.otherEvent.point;
      if (p1[0] === p2[0] && p1[1] === p2[1]) return 0;
      else return le1.contourId > le2.contourId ? 1 : -1;
    }
  } else { // Segments are collinear, but belong to separate polygons
    return le1.isSubject ? -1 : 1;
  }

  return compareEvents(le1, le2) === 1 ? 1 : -1;
}

function subdivide(eventQueue, subject, clipping, sbbox, cbbox, operation) {
  const sweepLine = new SplayTree(compareSegments);
  const sortedEvents = [];

  const rightbound = Math.min(sbbox[2], cbbox[2]);

  let prev, next, begin;

  while (eventQueue.length !== 0) {
    let event = eventQueue.pop();
    sortedEvents.push(event);

    // optimization by bboxes for intersection and difference goes here
    if ((operation === INTERSECTION && event.point[0] > rightbound) ||
        (operation === DIFFERENCE   && event.point[0] > sbbox[2])) {
      break;
    }

    if (event.left) {
      next  = prev = sweepLine.insert(event);
      begin = sweepLine.minNode();

      if (prev !== begin) prev = sweepLine.prev(prev);
      else                prev = null;

      next = sweepLine.next(next);

      const prevEvent = prev ? prev.key : null;
      let prevprevEvent;
      computeFields(event, prevEvent, operation);
      if (next) {
        if (possibleIntersection(event, next.key, eventQueue) === 2) {
          computeFields(event, prevEvent, operation);
          computeFields(event, next.key, operation);
        }
      }

      if (prev) {
        if (possibleIntersection(prev.key, event, eventQueue) === 2) {
          let prevprev = prev;
          if (prevprev !== begin) prevprev = sweepLine.prev(prevprev);
          else                    prevprev = null;

          prevprevEvent = prevprev ? prevprev.key : null;
          computeFields(prevEvent, prevprevEvent, operation);
          computeFields(event,     prevEvent,     operation);
        }
      }
    } else {
      event = event.otherEvent;
      next = prev = sweepLine.find(event);

      if (prev && next) {

        if (prev !== begin) prev = sweepLine.prev(prev);
        else                prev = null;

        next = sweepLine.next(next);
        sweepLine.remove(event);

        if (next && prev) {
          possibleIntersection(prev.key, next.key, eventQueue);
        }
      }
    }
  }
  return sortedEvents;
}

/**
 * @param  {Array.<SweepEvent>} sortedEvents
 * @return {Array.<SweepEvent>}
 */
function orderEvents(sortedEvents) {
  let event, i, len, tmp;
  const resultEvents = [];
  for (i = 0, len = sortedEvents.length; i < len; i++) {
    event = sortedEvents[i];
    if ((event.left && event.inResult) ||
      (!event.left && event.otherEvent.inResult)) {
      resultEvents.push(event);
    }
  }
  // Due to overlapping edges the resultEvents array can be not wholly sorted
  let sorted = false;
  while (!sorted) {
    sorted = true;
    for (i = 0, len = resultEvents.length; i < len; i++) {
      if ((i + 1) < len &&
        compareEvents(resultEvents[i], resultEvents[i + 1]) === 1) {
        tmp = resultEvents[i];
        resultEvents[i] = resultEvents[i + 1];
        resultEvents[i + 1] = tmp;
        sorted = false;
      }
    }
  }


  for (i = 0, len = resultEvents.length; i < len; i++) {
    event = resultEvents[i];
    event.pos = i;
  }

  // imagine, the right event is found in the beginning of the queue,
  // when his left counterpart is not marked yet
  for (i = 0, len = resultEvents.length; i < len; i++) {
    event = resultEvents[i];
    if (!event.left) {
      tmp = event.pos;
      event.pos = event.otherEvent.pos;
      event.otherEvent.pos = tmp;
    }
  }

  return resultEvents;
}


/**
 * @param  {Number} pos
 * @param  {Array.<SweepEvent>} resultEvents
 * @param  {Object>}    processed
 * @return {Number}
 */
function nextPos(pos, resultEvents, processed, origIndex) {
  let newPos = pos + 1;
  const length = resultEvents.length;
  if (newPos > length - 1) return pos - 1;
  let p  = resultEvents[pos].point;
  let p1 = resultEvents[newPos].point;


  // while in range and not the current one by value
  while (newPos < length && p1[0] === p[0] && p1[1] === p[1]) {
    if (!processed[newPos]) {
      return newPos;
    } else   {
      newPos++;
    }
    p1 = resultEvents[newPos].point;
  }

  newPos = pos - 1;

  while (processed[newPos] && newPos >= origIndex) {
    newPos--;
  }
  return newPos;
}


/**
 * @param  {Array.<SweepEvent>} sortedEvents
 * @return {Array.<*>} polygons
 */
function connectEdges(sortedEvents, operation) {
  let i, len;
  const resultEvents = orderEvents(sortedEvents);

  // "false"-filled array
  const processed = {};
  const result = [];
  let event;

  for (i = 0, len = resultEvents.length; i < len; i++) {
    if (processed[i]) continue;
    const contour = [[]];

    if (!resultEvents[i].isExteriorRing) {
      if (operation === DIFFERENCE && !resultEvents[i].isSubject && result.length === 0) {
        result.push(contour);
      } else if (result.length === 0) {
        result.push([[contour]]);
      } else {
        result[result.length - 1].push(contour[0]);
      }
    } else if (operation === DIFFERENCE && !resultEvents[i].isSubject && result.length > 1) {
      result[result.length - 1].push(contour[0]);
    } else {
      result.push(contour);
    }

    const ringId = result.length - 1;
    let pos = i;

    const initial = resultEvents[i].point;
    contour[0].push(initial);

    while (pos >= i) {
      event = resultEvents[pos];
      processed[pos] = true;

      if (event.left) {
        event.resultInOut = false;
        event.contourId   = ringId;
      } else {
        event.otherEvent.resultInOut = true;
        event.otherEvent.contourId   = ringId;
      }

      pos = event.pos;
      processed[pos] = true;
      contour[0].push(resultEvents[pos].point);
      pos = nextPos(pos, resultEvents, processed, i);
    }

    pos = pos === -1 ? i : pos;

    event = resultEvents[pos];
    processed[pos] = processed[event.pos] = true;
    event.otherEvent.resultInOut = true;
    event.otherEvent.contourId   = ringId;
  }

  // Handle if the result is a polygon (eg not multipoly)
  // Commented it again, let's see what do we mean by that
  // if (result.length === 1) result = result[0];
  return result;
}

var tinyqueue = TinyQueue;
var default_1 = TinyQueue;

function TinyQueue(data, compare) {
    if (!(this instanceof TinyQueue)) return new TinyQueue(data, compare);

    this.data = data || [];
    this.length = this.data.length;
    this.compare = compare || defaultCompare;

    if (this.length > 0) {
        for (var i = (this.length >> 1) - 1; i >= 0; i--) this._down(i);
    }
}

function defaultCompare(a, b) {
    return a < b ? -1 : a > b ? 1 : 0;
}

TinyQueue.prototype = {

    push: function (item) {
        this.data.push(item);
        this.length++;
        this._up(this.length - 1);
    },

    pop: function () {
        if (this.length === 0) return undefined;

        var top = this.data[0];
        this.length--;

        if (this.length > 0) {
            this.data[0] = this.data[this.length];
            this._down(0);
        }
        this.data.pop();

        return top;
    },

    peek: function () {
        return this.data[0];
    },

    _up: function (pos) {
        var data = this.data;
        var compare = this.compare;
        var item = data[pos];

        while (pos > 0) {
            var parent = (pos - 1) >> 1;
            var current = data[parent];
            if (compare(item, current) >= 0) break;
            data[pos] = current;
            pos = parent;
        }

        data[pos] = item;
    },

    _down: function (pos) {
        var data = this.data;
        var compare = this.compare;
        var halfLength = this.length >> 1;
        var item = data[pos];

        while (pos < halfLength) {
            var left = (pos << 1) + 1;
            var right = left + 1;
            var best = data[left];

            if (right < this.length && compare(data[right], best) < 0) {
                left = right;
                best = data[right];
            }
            if (compare(best, item) >= 0) break;

            data[pos] = best;
            pos = left;
        }

        data[pos] = item;
    }
};
tinyqueue.default = default_1;

const max = Math.max;
const min = Math.min;

let contourId = 0;


function processPolygon(contourOrHole, isSubject, depth, Q, bbox, isExteriorRing) {
  let i, len, s1, s2, e1, e2;
  for (i = 0, len = contourOrHole.length - 1; i < len; i++) {
    s1 = contourOrHole[i];
    s2 = contourOrHole[i + 1];
    e1 = new SweepEvent(s1, false, undefined, isSubject);
    e2 = new SweepEvent(s2, false, e1,        isSubject);
    e1.otherEvent = e2;

    if (s1[0] === s2[0] && s1[1] === s2[1]) {
      continue; // skip collapsed edges, or it breaks
    }

    e1.contourId = e2.contourId = depth;
    if (!isExteriorRing) {
      e1.isExteriorRing = false;
      e2.isExteriorRing = false;
    }
    if (compareEvents(e1, e2) > 0) {
      e2.left = true;
    } else {
      e1.left = true;
    }

    const x = s1[0], y = s1[1];
    bbox[0] = min(bbox[0], x);
    bbox[1] = min(bbox[1], y);
    bbox[2] = max(bbox[2], x);
    bbox[3] = max(bbox[3], y);

    // Pushing it so the queue is sorted from left to right,
    // with object on the left having the highest priority.
    Q.push(e1);
    Q.push(e2);
  }
}


function fillQueue(subject, clipping, sbbox, cbbox, operation) {
  const eventQueue = new tinyqueue(null, compareEvents);
  let polygonSet, isExteriorRing, i, ii, j, jj; //, k, kk;

  for (i = 0, ii = subject.length; i < ii; i++) {
    polygonSet = subject[i];
    for (j = 0, jj = polygonSet.length; j < jj; j++) {
      isExteriorRing = j === 0;
      if (isExteriorRing) contourId++;
      processPolygon(polygonSet[j], true, contourId, eventQueue, sbbox, isExteriorRing);
    }
  }

  for (i = 0, ii = clipping.length; i < ii; i++) {
    polygonSet = clipping[i];
    for (j = 0, jj = polygonSet.length; j < jj; j++) {
      isExteriorRing = j === 0;
      if (operation === DIFFERENCE) isExteriorRing = false;
      if (isExteriorRing) contourId++;
      processPolygon(polygonSet[j], false, contourId, eventQueue, cbbox, isExteriorRing);
    }
  }

  return eventQueue;
}

const EMPTY = [];


function trivialOperation(subject, clipping, operation) {
  let result = null;
  if (subject.length * clipping.length === 0) {
    if        (operation === INTERSECTION) {
      result = EMPTY;
    } else if (operation === DIFFERENCE) {
      result = subject;
    } else if (operation === UNION ||
               operation === XOR) {
      result = (subject.length === 0) ? clipping : subject;
    }
  }
  return result;
}


function compareBBoxes(subject, clipping, sbbox, cbbox, operation) {
  let result = null;
  if (sbbox[0] > cbbox[2] ||
      cbbox[0] > sbbox[2] ||
      sbbox[1] > cbbox[3] ||
      cbbox[1] > sbbox[3]) {
    if        (operation === INTERSECTION) {
      result = EMPTY;
    } else if (operation === DIFFERENCE) {
      result = subject;
    } else if (operation === UNION ||
               operation === XOR) {
      result = subject.concat(clipping);
    }
  }
  return result;
}


function boolean(subject, clipping, operation) {
  if (typeof subject[0][0][0] === 'number') {
    subject = [subject];
  }
  if (typeof clipping[0][0][0] === 'number') {
    clipping = [clipping];
  }
  let trivial = trivialOperation(subject, clipping, operation);
  if (trivial) {
    return trivial === EMPTY ? null : trivial;
  }
  const sbbox = [Infinity, Infinity, -Infinity, -Infinity];
  const cbbox = [Infinity, Infinity, -Infinity, -Infinity];

  //console.time('fill queue');
  const eventQueue = fillQueue(subject, clipping, sbbox, cbbox, operation);
  //console.timeEnd('fill queue');

  trivial = compareBBoxes(subject, clipping, sbbox, cbbox, operation);
  if (trivial) {
    return trivial === EMPTY ? null : trivial;
  }
  //console.time('subdivide edges');
  const sortedEvents = subdivide(eventQueue, subject, clipping, sbbox, cbbox, operation);
  //console.timeEnd('subdivide edges');

  //console.time('connect vertices');
  const result = connectEdges(sortedEvents, operation);
  //console.timeEnd('connect vertices');
  return result;
}

var intersect_1 = createCommonjsModule(function (module, exports) {
var __importStar = (commonjsGlobal && commonjsGlobal.__importStar) || function (mod) {
    if (mod && mod.__esModule) return mod;
    var result = {};
    if (mod != null) for (var k in mod) if (Object.hasOwnProperty.call(mod, k)) result[k] = mod[k];
    result["default"] = mod;
    return result;
};
Object.defineProperty(exports, "__esModule", { value: true });


var martinez = __importStar(boolean);
/**
 * Takes two {@link Polygon|polygon} or {@link MultiPolygon|multi-polygon} geometries and
 * finds their polygonal intersection. If they don't intersect, returns null.
 *
 * @name intersect
 * @param {Feature<Polygon | MultiPolygon>} poly1 the first polygon or multipolygon
 * @param {Feature<Polygon | MultiPolygon>} poly2 the second polygon or multipolygon
 * @param {Object} [options={}] Optional Parameters
 * @param {Object} [options.properties={}] Translate GeoJSON Properties to Feature
 * @returns {Feature|null} returns a feature representing the area they share (either a {@link Polygon} or
 * {@link MultiPolygon}). If they do not share any area, returns `null`.
 * @example
 * var poly1 = turf.polygon([[
 *   [-122.801742, 45.48565],
 *   [-122.801742, 45.60491],
 *   [-122.584762, 45.60491],
 *   [-122.584762, 45.48565],
 *   [-122.801742, 45.48565]
 * ]]);
 *
 * var poly2 = turf.polygon([[
 *   [-122.520217, 45.535693],
 *   [-122.64038, 45.553967],
 *   [-122.720031, 45.526554],
 *   [-122.669906, 45.507309],
 *   [-122.723464, 45.446643],
 *   [-122.532577, 45.408574],
 *   [-122.487258, 45.477466],
 *   [-122.520217, 45.535693]
 * ]]);
 *
 * var intersection = turf.intersect(poly1, poly2);
 *
 * //addToMap
 * var addToMap = [poly1, poly2, intersection];
 */
function intersect(poly1, poly2, options) {
    if (options === void 0) { options = {}; }
    var geom1 = invariant.getGeom(poly1);
    var geom2 = invariant.getGeom(poly2);
    if (geom1.type === "Polygon" && geom2.type === "Polygon") {
        var intersection = martinez.intersection(geom1.coordinates, geom2.coordinates);
        if (intersection === null || intersection.length === 0) {
            return null;
        }
        if (intersection.length === 1) {
            var start = intersection[0][0][0];
            var end = intersection[0][0][intersection[0][0].length - 1];
            if (start[0] === end[0] && start[1] === end[1]) {
                return helpers.polygon(intersection[0], options.properties);
            }
            return null;
        }
        return helpers.multiPolygon(intersection, options.properties);
    }
    else if (geom1.type === "MultiPolygon") {
        var resultCoords = [];
        // iterate through the polygon and run intersect with each part, adding to the resultCoords.
        for (var _i = 0, _a = geom1.coordinates; _i < _a.length; _i++) {
            var coords = _a[_i];
            var subGeom = invariant.getGeom(helpers.polygon(coords));
            var subIntersection = intersect(subGeom, geom2);
            if (subIntersection) {
                var subIntGeom = invariant.getGeom(subIntersection);
                if (subIntGeom.type === "Polygon") {
                    resultCoords.push(subIntGeom.coordinates);
                }
                else if (subIntGeom.type === "MultiPolygon") {
                    resultCoords = resultCoords.concat(subIntGeom.coordinates);
                }
                else {
                    throw new Error("intersection is invalid");
                }
            }
        }
        // Make a polygon with the result
        if (resultCoords.length === 0) {
            return null;
        }
        if (resultCoords.length === 1) {
            return helpers.polygon(resultCoords[0], options.properties);
        }
        else {
            return helpers.multiPolygon(resultCoords, options.properties);
        }
    }
    else if (geom2.type === "MultiPolygon") {
        // geom1 is a polygon and geom2 a multiPolygon,
        // put the multiPolygon first and fallback to the previous case.
        return intersect(geom2, geom1);
    }
    else {
        // handle invalid geometry types
        throw new Error("poly1 and poly2 must be either polygons or multiPolygons");
    }
}
exports.default = intersect;
});

var intersect = unwrapExports(intersect_1);

function _getKrigingGridInfo(featureCollection,weight,krigingParams){
    //先获取featureCollection的bbox
    let values=[],lons=[],lats=[];
    let extent=[100000000,100000000,-100000000,-100000000];
    featureCollection.features.forEach(feature => {
        //提取插值权重字段，准备克里金插值使用
        values.push(feature.properties[weight]);
        lons.push(feature.geometry.coordinates[0]);
        lats.push(feature.geometry.coordinates[1]);
        if(extent[0]>feature.geometry.coordinates[0])
            extent[0]=feature.geometry.coordinates[0];
        if(extent[1]>feature.geometry.coordinates[1])
            extent[1]=feature.geometry.coordinates[1];
        if(extent[2]<feature.geometry.coordinates[0])
            extent[2]=feature.geometry.coordinates[0];
        if(extent[3]<feature.geometry.coordinates[1])
            extent[3]=feature.geometry.coordinates[1];
    });
    let variogram=kriging.train(values,lons,lats,krigingParams.model,krigingParams.sigma2,krigingParams.alpha);
    let gridinfo=kriging.getGridInfo(extent,variogram,200);
    return gridinfo;
}

/*
* 克里金生成矢量等值面，浏览器和node都可以使用
* @param {json} featureCollection：必填，已有点数据，geojson格式
* @param {string} weight：必填，插值所依赖的圈中字段
* @param {object) krigingParams：必填，克里金插值算法参数设置
    krigingParams:{
         krigingModel:'exponential','gaussian','spherical'，三选一
         krigingSigma2:
         krigingAlpha:
    }
* @param {array} breaks：必填，等值面分级区间
* @param {json) clip_feature：可选，切割面要素的geojson格式
*/
function getVectorContour(featureCollection,weight,krigingParams,breaks,clip_feature){
    let gridinfo=_getKrigingGridInfo(featureCollection,weight,krigingParams);
    let vectorContour=kriging.getVectorContour(gridinfo,breaks);
     //是否需要切割
     if(clip_feature){
		let clip_features=[];
        vectorContour.features.forEach(feature=>{
			let _feature=intersect(feature,clip_feature);
			//补全切割要素属性信息
			if(_feature){
				_feature.properties=feature.properties;
				clip_features.push(_feature);
			}
        });
        vectorContour.features=clip_features;
    }
    return vectorContour;
}
/*
* 克里金生成栅格等值面并绘制到canvas上，仅浏览器中使用
* @param {json} featureCollection：必填，已有点数据，geojson格式
* @param {string} weight：必填，插值所依赖的圈中字段
* @param {object) krigingParams：必填，克里金插值算法参数设置
    krigingParams:{
         krigingModel:'exponential','gaussian','spherical'，三选一
         krigingSigma2:
         krigingAlpha:
    }
* @param {dom) canvas：必填，绑定渲染的canvas dom
* @param {array) colors：必填，等值面分级区间
*/
function drawCanvasContour(featureCollection,weight,krigingParams,canvas,xlim,ylim,colors) {
	let gridinfo=_getKrigingGridInfo(featureCollection,weight,krigingParams);
    kriging.drawCanvasContour(gridinfo,canvas,xlim,ylim,colors);
}

exports.drawCanvasContour = drawCanvasContour;
exports.getVectorContour = getVectorContour;

Object.defineProperty(exports, '__esModule', { value: true });

})));

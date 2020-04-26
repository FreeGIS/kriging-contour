// kriging-contour v1.0.4 Copyright 2020 freegis
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
*/
function getVectorContour(featureCollection,weight,krigingParams,breaks){
    let gridinfo=_getKrigingGridInfo(featureCollection,weight,krigingParams);
    let vectorContour=kriging.getVectorContour(gridinfo,breaks);
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

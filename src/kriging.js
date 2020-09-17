import {contours as d3_contours} from 'd3-contour';

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
}

// Matrix algebra
function kriging_matrix_diag(c, n) {
	var i,
	Z = [0].rep(n * n);
	for (i = 0; i < n; i++)
		Z[i * n + i] = c;
	return Z;
};
function kriging_matrix_transpose(X, n, m) {
	var i,
	j,
	Z = Array(m * n);
	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
			Z[j * n + i] = X[i * m + j];
	return Z;
};
function kriging_matrix_scale(X, c, n, m) {
	var i,
	j;
	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
			X[i * m + j] *= c;
};
function kriging_matrix_add(X, Y, n, m) {
	var i,
	j,
	Z = Array(n * m);
	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
			Z[i * m + j] = X[i * m + j] + Y[i * m + j];
	return Z;
};
// Naive matrix multiplication
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
};
// Cholesky decomposition
function kriging_matrix_chol(X, n) {
	var i,
	j,
	k,
	sum,
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
};
// Inversion of cholesky decomposition
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

};
// Inversion via gauss-jordan elimination
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
};
function kriging_variogram_exponential(h, nugget, range, sill, A) {
	return nugget + ((sill - nugget) / range) *
	(1.0 - Math.exp( - (1.0 / A) * (h / range)));
};
function kriging_variogram_spherical(h, nugget, range, sill, A) {
	if (h > range)
		return nugget + (sill - nugget) / range;
	return nugget + ((sill - nugget) / range) *
	(1.5 * (h / range) - 0.5 * Math.pow(h / range, 3));
};

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
	};

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
		};
		Y[i] = semi[i];
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
	for(let i=0;i<_contours.length;i++){
		const contour = _contours[i];
		if(contour.type==='MultiPolygon'){
			contour.coordinates.forEach(polygon => {
				let geom = {
					"type" : "Polygon",
					"coordinates" : []
				};
				//坐标转换，图形为空去除
				geom.coordinates=polygon_pixel2geos(polygon,gridInfo);
				if(geom.coordinates.length>0){
					dataset.features.push({
						"type" : "Feature",
						"properties" : {
							"value" : contour.value
						},
						"geometry" :geom
					});
				}
			});
		}
		else if(contour.type==='Polygon'){
			let geom = {
				"type" : "Polygon",
				"coordinates" : []
			};
			//坐标转换，图形为空去除
			geom.coordinates=polygon_pixel2geos(contour.coordinates,gridInfo);
			if(geom.coordinates.length>0){
				dataset.features.push({
					"type" : "Feature",
					"properties" : {
						"value" : contour.value
					},
					"geometry" :geom
				});
			}
		}
	}
	return dataset;
};
//像素坐标转地理坐标
function polygon_pixel2geos(polygon,gridInfo){
	//polygon分内环和外环
	const _polygon = polygon.map((ring) => {
		const _ring = ring.map(function (coor) {
			//像素坐标转地理坐标 ，像素坐标y方向从上到下递增，纬度是y从上到下递减
			const lon = gridInfo.xlim[0] + coor[0]*gridInfo.x_resolution;
			let lat;
			//格网自上向下走
			if(gridInfo.y_resolution<0)
				lat=gridInfo.ylim[1] + coor[1]*gridInfo.y_resolution;
			//格网自下向上走
			else
				lat=gridInfo.ylim[0] + coor[1]*gridInfo.y_resolution;
	
			return [lon,lat];
		});
		return _ring;
	});	
	return _polygon;
}
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

export {kriging};

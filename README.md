# kriging-contour
基于克里金插值算法，根据离散点位置及其权重，生成等值面矢量数据(GeoJSON格式)和栅格数据(Canvas绘制图片)，这些数据在任何WebGIS客户端上都可通用展示。

## 兼容性

该库生成结果是geojson和canvas上绘制两种形式，这种通用结果可以在ol,mapbox,leaflet,arcgis api等各种前端api上直接展示，具有跨库的通用性。

## 安装
### 源码下载并编译：
```
# 下载源码安装依赖
git clone git@github.com:FreeGIS/kriging-contour.git
cd kriging-contour
npm install
# 编译
npm run build
```
### npm安装
```
npm install kriging-contour --save
```
## 使用说明
​	**矢量等值面**

kriging-contour.getVectorContour(dataset,weight_field,kriging_params,weight_breaks,clip_feature);

```
dataset：geojson格式的featureclass数据集，feature是图形是点，必填
weight_field：绑定权重字段名称，必填
kriging_params：克里金插值参数，必填
weight_breaks：权重生成等值面分级数组，必填
clip_feature: 切割范围，geojson格式的面要素，可选参数
```
​		示例代码：
```
	//计算克里金等值面
		let kriging_contours=kriging-contour.getVectorContourr(dataset,'level',{
			model:'exponential',
			sigma2:0,
			alpha:100
		},[0,10,20,30,40,50,60,70,80,90,100]);
```
​	**图片等值面**

kriging-contour.drawCanvasContour(dataset,weight_field,kriging_params,weight_breaks,clip_feature);
```
dataset：geojson格式的featureclass数据集，feature是图形是点，必填
weight_field：绑定权重字段名称，必填
kriging_params：克里金插值参数，必填
canvas：渲染的canvas对象，必填
xlim：当前视图窗口(extent)的x轴跨度。
ylim：当前视图窗口(extent)的y轴跨度。
colors：渲染颜色分级。
```
示例代码：
```
	//计算克里金等值面
		kriging-contour.drawCanvasContour(dataset,'level',{
			model:'exponential',
			sigma2:0,
			alpha:100
		},canvas,[extent[0],extent[2]],[extent[1],extent[3]],params.colors);

```



kriging是基于oeo4b的[kriging.js](https://github.com/oeo4b/kriging.js)修改的，原kriging.js编码不够规范，且仅仅支持将插值结果渲染到canvas形成等值面，这种图片形式的等值面锯齿比较严重。

本次修改修复部分编码不可读部分，且重新实现生成基于矢量的插值等值面，渲染效果较好。

kriging图片渲染效果:
![kriging图片渲染效果](https://github.com/FreeGIS/kriging-contour/blob/master/doc/raster.jpg)
kriging矢量渲染效果:
![kriging矢量渲染效果](https://github.com/FreeGIS/kriging-contour/blob/master/doc/vector.jpg)
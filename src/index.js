import {kriging} from './kriging';
import intersect from '@turf/intersect';

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
};

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
};
export {getVectorContour,drawCanvasContour};
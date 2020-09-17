import commonjs from '@rollup/plugin-commonjs';
import json from '@rollup/plugin-json';
import resolve from '@rollup/plugin-node-resolve';
import babel from 'rollup-plugin-babel';
import sourcemaps from 'rollup-plugin-sourcemaps';
import {terser} from 'rollup-plugin-terser';

export default {
    input: 'src/index.js',
    output: {
      name: 'kriging-contour',
      file: 'dist/kriging-contour.js',
      format: 'umd',
      sourcemap: true
    },
    treeshake: true,
    plugins: [ 
      json(),
      resolve(),
      commonjs(),
      sourcemaps(),
      babel({
        exclude: 'node_modules/**' // 只编译我们的源代码
      }),
      terser()
    ]
  };
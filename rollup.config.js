import * as meta from "./package.json";
import commonjs from 'rollup-plugin-commonjs';
import resolve from 'rollup-plugin-node-resolve';
import {terser} from "rollup-plugin-terser";
const config = {
  input: "src/index.js",
  output: {
    file: `dist/${meta.name}.js`,
    name: "kriging",
    format: "umd",
    indent: false,
    extend: true,
    banner: `// ${meta.name} v${meta.version} Copyright ${(new Date).getFullYear()} ${meta.author}`
  },
  plugins: [resolve(),commonjs()]
};
export default [
  config,
  {
    ...config,
    output: {
      ...config.output,
      file: `dist/${meta.name}.min.js`
    },
    plugins: [
      ...config.plugins,
      terser({
        output: {
          preamble: config.output.banner
        }
      })
    ]
  }
];
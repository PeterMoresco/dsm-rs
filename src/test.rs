// extern crate ndarray;
// use ndarray::prelude::*;

fn ff(a: f64) {
    println!("This is {}", a);
}

fn main() {
//     //let a = [1, 2, 3, 4, 5];
    // let nos: [[f64; 2]; 6] = [[0.0,0.0], [0.0, 9.144], [9.144, 0.0],[9.144,9.144],[18.288,0.0],[18.288,9.144]];
    // let nos: Vec<[f64;2]> = vec![[0.0,0.0], [0.0, 9.144], [9.144, 0.0],[9.144,9.144],[18.288,0.0],[18.288,9.144]];
    // let barras: Vec<[i32; 2]> = vec![[0,2],[0,3],[1,2],[1,3],[3,2],[3,4],[5,2],[3,5],[2,4],[4,5]];
//     let temp: f64 = &nos[3].unwrap().clone();
//     // temp.what_is_this();
//     let b: f64 = 5.0;
//     println!("{}", b.sqrt());
//     //println!("{:?}", nos[1][0]);
//     //barras.what_is_this();
//     let mut v = Vec::new();
//     let mut _f = 0;

//     for _ in nos.iter() {
//         v.push((_f, _f + 1));
//         _f =  _f + 2;
//     }

//     for (ind, bar) in barras.iter().enumerate() {
//         println!("{}", ind);
//         println!("{:?}", bar);

//     }
let a = 5.;
ff(a);
println!("This is from outside the fn {}", a);
}
/*
This is a re-implementation of the DS_Method.py
originally made by me for my bachelor thesis.
It's basically a calculator for the direct stiffnes method for
2D trusses.
*/
#![allow(non_snake_case)]
use std::vec::Vec;
// use ndarray::{ Array3, Array, Array2, arr2, array };
// use ndarray:: prelude::*;
use ndarray::*;
use ndarray_linalg::*;


fn stiffness_matrix(nos: &Vec<[f64; 2]>, barras: &Vec<[i32; 2]>, areas: &Vec<f64>, mom_e: f64) -> (Array2<f64>, Array3<f64>) {
    // Here the 2D vectors are associated with 
    // the nodes
    let mut f_index = Vec::new();
    let mut _f = 0;
    // Populate the matrix
    //TODO this could be improved
    for _ in nos.iter() {
        f_index.push([_f, _f + 1]);
        _f = _f + 2
    }
    // Create the matrix to hold all the 
    // local bars matrix
    let mut mrl: Array3<f64> = Array::zeros((4, 4, barras.len()));
    // Create the matrix to hold all bars lenghts
    let mut comp_barras: Vec<f64> = Vec::new();
    // Calculate bars length
    for i in 0..barras.len(){
        let _x1: f64 = nos[barras[i][0] as usize][0];
        let _x2: f64 = nos[barras[i][1] as usize][0];
        let _y1: f64 = nos[barras[i][0] as usize][1];
        let _y2: f64 = nos[barras[i][1] as usize][1];
        let _len: f64 = (_x2-_x1).hypot(_y2-_y1);
        comp_barras.push(_len);
    }
    // Calculate the freedom degree
    let gL: usize = 2 * nos.len();
    // Creates the global stiffness matrix
    let mut mrg: Array2<f64> = Array::zeros((gL, gL));
    // Initiate the loop trough the bars and
    // calculate the local stiffness matrices
    for (ind, bar) in barras.iter().enumerate() {
        // Converts bar nodes to usize
        let _n0: usize = bar[0] as usize;
        let _n1: usize = bar[1] as usize;
        // Calculates the cossene vectors
        let Cx: f64 = (nos[_n1][0] - nos[_n0][0])/comp_barras[ind];
        let Cy: f64 = (nos[_n1][1] - nos[_n0][1])/comp_barras[ind];
        // Temp matrix to ease index
        let mut mrlb: Array2<f64> = Array::zeros((4, 4));
        // Fillss the first quadrant
        mrlb[[0, 0]] = Cx * Cx;
        mrlb[[1, 1]] = Cy * Cy;
        mrlb[[1, 0]] = Cx * Cy;
        mrlb[[0, 1]] = Cx * Cy;
        // Fills the second quadrant
        let _tmp_slice: Array2<f64> = -1. * mrlb.slice(s![0..2, 0..2]).to_owned();
        mrlb.slice_mut(s![2..4, 0..2]).assign(&_tmp_slice);
        // Fills the third and fourth quadrant
        let _tmp_slice: Array2<f64> = -1. * mrlb.slice(s![.., 0..2]).to_owned();
        mrlb.slice_mut(s![.., 2..4]).assign(&_tmp_slice);
        // println!("This is the bar {} and \n {:?}", &ind+1, &mrlb);
        // Multiplies the values by the areas
        mrlb = mrlb * areas[ind];
        // Multiplies by the length of th bars
        mrlb = mrlb / comp_barras[ind];
        // println!("This is the {} bar area {} and l {}", ind+1, areas[ind], comp_barras[ind]);
        // Create the array of the forces by the nodes
        let forca = [f_index[_n0][0],
                     f_index[_n0][1],
                     f_index[_n1][0],
                     f_index[_n1][1]];
        // Keep the identity of the forces
        let mut mat_id: Array2<f64> = Array::zeros((gL, 4));
        for i in 0..4 {
            mat_id[[forca[i], i]] = 1.;
        }
        let prev = mat_id.dot(&mrlb);
        mrg = mrg + prev.dot(&mat_id.t());
        mrl.slice_mut(s![.., .., ind]).assign(&mrlb);
    };
    mrg = mrg * mom_e;
    mrl = mrl * mom_e;
    (mrg, mrl)
}
// Implement the solve methd https://notesformsc.org/c-gauss-elimination-method/

fn main() {
    // Momento de inércia
    let mom: f64 = 68950000.0;
    // Nós da treliça
    let nos: Vec<[f64; 2]> = vec![[0.0,0.0], [0.0, 9.144], [9.144, 0.0],[9.144,9.144],[18.288,0.0],[18.288,9.144]];
    // Barras da treliça
    let barras: Vec<[i32; 2]> = vec![[0,2],[0,3],[1,2],[1,3],[3,2],[3,4],[5,2],[3,5],[2,4],[4,5]];
    // Areas das secoes
    let areas: Vec<f64> = vec![19519.11, 64.5, 14784.81, 9705.02, 65.48,
                               300.48, 4819.16, 13610.37, 13930.3, 64.5];
    // Carregamento
    // let n_load = [5, 9];
    // let f_load = [-444.82, -444.82];
    // Nós com restrição
    // let drest = [0,1,2,3];
    let a = stiffness_matrix(&nos, &barras, &areas, mom);
    // println!("This is the mrg\n{:?}", a.0)
    let b: Array2<f64> = random((3, 3));
    let c: Array1<f64> = random(3);
    let _x = b.solve(&c);
}
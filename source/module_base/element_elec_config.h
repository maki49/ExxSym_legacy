#ifndef ELEMENT_ELEC_CONFIG
#define ELEMENT_ELEC_CONFIG

#include <map>
#include <string>

namespace ModuleBase
{

const std::map<std::string, std::string> EleConfig
= {
{"H",  "1s1"},   
{"He", "1s2"},  
{"Li", "[He] 2s1"},  
{"Be", "[He] 2s2"},  
{"B",  "[He] 2s2 2p1"},   
{"C",  "[He] 2s2 2p2"},   
{"N",  "[He] 2s2 2p3"},   
{"O",  "[He] 2s2 2p4"},   
{"F",  "[He] 2s2 2p5"},
{"Ne", "[He] 2s2 2p6"}, 
{"Na", "[Ne] 3s1"}, 
{"Mg", "[Ne] 3s2"}, 
{"Al", "[Ne] 3s2 3p1"}, 
{"Si", "[Ne] 3s2 3p2"}, 
{"P",  "[Ne] 3s2 3p3"},  
{"S",  "[Ne] 3s2 3p4"},  
{"Cl", "[Ne] 3s2 3p5"}, 
{"Ar", "[Ne] 3s2 3p6"},
{"K",  "[Ar] 4s1"},  
{"Ca", "[Ar] 4s2"}, 
{"Sc", "[Ar] 3d1 4s2"}, 
{"Ti", "[Ar] 3d2 4s2"}, 
{"V",  "[Ar] 3d3 4s2"},  
{"Cr", "[Ar] 3d5 4s1"}, 
{"Mn", "[Ar] 3d5 4s2"}, 
{"Fe", "[Ar] 3d6 4s2"}, 
{"Co", "[Ar] 3d7 4s2"},
{"Ni", "[Ar] 3d8 4s2"}, 
{"Cu", "[Ar] 3d10 4s1"}, 
{"Zn", "[Ar] 3d10 4s2"}, 
{"Ga", "[Ar] 3d10 4s2 4p1"}, 
{"Ge", "[Ar] 3d10 4s2 4p2"}, 
{"As", "[Ar] 3d10 4s2 4p3"}, 
{"Se", "[Ar] 3d10 4s2 4p4"}, 
{"Br", "[Ar] 3d10 4s2 4p5"}, 
{"Kr", "[Ar] 3d10 4s2 4p6"},
{"Rb", "[Kr] 5s1"}, 
{"Sr", "[Kr] 5s2"}, 
{"Y",  "[Kr] 4d1 5s2"}, 
{"Zr", "[Kr] 4d2 5s2"}, 
{"Nb", "[Kr] 4d4 5s1"}, 
{"Mo", "[Kr] 4d5 5s1"}, 
{"Tc", "[Kr] 4d5 5s2"}, 
{"Ru", "[Kr] 4d7 5s1"}, 
{"Rh", "[Kr] 4d8 5s1"},
{"Pd", "[Kr] 4d10"}, 
{"Ag", "[Kr] 4d10 5s1"}, 
{"Cd", "[Kr] 4d10 5s2"}, 
{"In", "[Kr] 4d10 5s2 5p1"}, 
{"Sn", "[Kr] 4d10 5s2 5p2"}, 
{"Sb", "[Kr] 4d10 5s2 5p3"}, 
{"Te", "[Kr] 4d10 5s2 5p4"}, 
{"I",  "[Kr] 4d10 5s2 5p5"},  
{"Xe", "[Kr] 4d10 5s2 5p6"},
{"Cs", "[Xe] 6s1"}, 
{"Ba", "[Xe] 6s2"}, 
{"La", "[Xe] 5d1 6s2"}, 
{"Ce", "[Xe] 4f1 5d1 6s2"}, 
{"Pr", "[Xe] 4f3 6s2"}, 
{"Nd", "[Xe] 4f4 6s2"}, 
{"Pm", "[Xe] 4f5 6s2"}, 
{"Sm", "[Xe] 4f6 6s2"}, 
{"Eu", "[Xe] 4f7 6s2"},
{"Gd", "[Xe] 4f7 5d1 6s2"}, 
{"Tb", "[Xe] 4f9 6s2"}, 
{"Dy", "[Xe] 4f10 6s2"}, 
{"Ho", "[Xe] 4f11 6s2"}, 
{"Er", "[Xe] 4f12 6s2"}, 
{"Tm", "[Xe] 4f13 6s2"}, 
{"Yb", "[Xe] 4f14 6s2"}, 
{"Lu", "[Xe] 4f14 5d1 6s2"}, 
{"Hf", "[Xe] 4f14 5d2 6s2"},
{"Ta", "[Xe] 4f14 5d3 6s2"}, 
{"W",  "[Xe] 4f14 5d4 6s2"},  
{"Re", "[Xe] 4f14 5d5 6s2"}, 
{"Os", "[Xe] 4f14 5d6 6s2"}, 
{"Ir", "[Xe] 4f14 5d7 6s2"}, 
{"Pt", "[Xe] 4f14 5d9 6s1"}, 
{"Au", "[Xe] 4f14 5d10 6s1"}, 
{"Hg", "[Xe] 4f14 5d10 6s2"}, 
{"Tl", "[Xe] 4f14 5d10 6s2 6p1"},
{"Pb", "[Xe] 4f14 5d10 6s2 6p2"}, 
{"Bi", "[Xe] 4f14 5d10 6s2 6p3"}, 
{"Po", "[Xe] 4f14 5d10 6s2 6p4"}, 
{"At", "[Xe] 4f14 5d10 6s2 6p5"}, 
{"Rn", "[Xe] 4f14 5d10 6s2 6p6"}, 
{"Fr", "[Rn] 7s1"},
{"Ra", "[Rn] 7s2"},
{"Ac", "[Rn] 6d1 7s2"},
{"Th", "[Rn] 6d2 7s2"},
{"Pa", "[Rn] 5f2 6d1 7s2"},
{"U" , "[Rn] 5f3 6d1 7s2"},
{"Np", "[Rn] 5f4 6d1 7s2"},
{"Pu", "[Rn] 5f6 7s2"},
{"Am", "[Rn] 5f7 7s2"},
{"Cm", "[Rn] 5f7 6d1 7s2"},
{"Bk", "[Rn] 5f9 7s2"},
{"Cf", "[Rn] 5f10 7s2"},
{"Es", "[Rn] 5f11 7s2"},
{"Fm", "[Rn] 5f12 7s2"},
{"Md", "[Rn] 5f13 7s2"},
{"No", "[Rn] 5f14 7s2"},
{"Lr", "[Rn] 5f14 7s2 7p1"},
{"Rf", "[Rn] 5f14 6d2 7s2"},
{"Db", "[Rn] 5f14 6d3 7s2"},
{"Sg", "[Rn] 5f14 6d4 7s2"},
{"Bh", "[Rn] 5f14 6d5 7s2"},
{"Hs", "[Rn] 5f14 6d6 7s2"},
{"Mt", "[Rn] 5f14 6d7 7s2"},
{"Ds", "[Rn] 5f14 6d8 7s2"},
{"Rg", "[Rn] 5f14 6d10 7s1"},
{"Cn", "[Rn] 5f14 6d10 7s2"},
{"Nh", "[Rn] 5f14 6d10 7s2 7p1"},
{"Fl", "[Rn] 5f14 6d10 7s2 7p2"},
{"Mc", "[Rn] 5f14 6d10 7s2 7p3"},
{"Lv", "[Rn] 5f14 6d10 7s2 7p4"},
{"Ts", "[Rn] 5f14 6d10 7s2 7p5"},
{"Og", "[Rn] 5f14 6d10 7s2 7p6"}
};

const std::map<std::string, int> MinZval
= {
{"H",  1},   
{"He", 2},  
{"Li", 1},  
{"Be", 2},  
{"B",  3},   
{"C",  4},   
{"N",  5},   
{"O",  6},   
{"F",  7},
{"Ne", 8}, 
{"Na", 1}, 
{"Mg", 2}, 
{"Al", 3}, 
{"Si", 4}, 
{"P",  5},  
{"S",  6},  
{"Cl", 7}, 
{"Ar", 8},
{"K",  1},  
{"Ca", 2}, 
{"Sc", 3}, 
{"Ti", 4}, 
{"V",  5},  
{"Cr", 6}, 
{"Mn", 7}, 
{"Fe", 8}, 
{"Co", 9},
{"Ni", 10}, 
{"Cu", 11}, 
{"Zn", 12}, 
{"Ga", 3}, 
{"Ge", 4}, 
{"As", 5}, 
{"Se", 6}, 
{"Br", 7}, 
{"Kr", 8},
{"Rb", 1}, 
{"Sr", 2}, 
{"Y",  3}, 
{"Zr", 4}, 
{"Nb", 5}, 
{"Mo", 6}, 
{"Tc", 7}, 
{"Ru", 8}, 
{"Rh", 9},
{"Pd", 10}, 
{"Ag", 11}, 
{"Cd", 12}, 
{"In", 3}, 
{"Sn", 4}, 
{"Sb", 5}, 
{"Te", 6}, 
{"I",  7},  
{"Xe", 8},
{"Cs", 1}, 
{"Ba", 2}, 
{"La", 3}, 
{"Ce", 4}, 
{"Pr", 5}, 
{"Nd", 6}, 
{"Pm", 7}, 
{"Sm", 8}, 
{"Eu", 9},
{"Gd", 10}, 
{"Tb", 11}, 
{"Dy", 12}, 
{"Ho", 13}, 
{"Er", 14}, 
{"Tm", 15}, 
{"Yb", 16}, 
{"Lu", 17}, 
{"Hf", 4},
{"Ta", 5}, 
{"W",  6},  
{"Re", 7}, 
{"Os", 8}, 
{"Ir", 9}, 
{"Pt", 10}, 
{"Au", 11}, 
{"Hg", 12}, 
{"Tl", 3},
{"Pb", 4}, 
{"Bi", 5}, 
{"Po", 6}, 
{"At", 7}, 
{"Rn", 8}, 
{"Fr", 1},
{"Ra", 2},
{"Ac", 3},
{"Th", 4},
{"Pa", 5},
{"U" , 6},
{"Np", 7},
{"Pu", 8},
{"Am", 9},
{"Cm", 10},
{"Bk", 11},
{"Cf", 12},
{"Es", 13},
{"Fm", 14},
{"Md", 15},
{"No", 16},
{"Lr", 17},
{"Rf", 18},
{"Db", 19},
{"Sg", 20},
{"Bh", 21},
{"Hs", 22},
{"Mt", 23},
{"Ds", 24},
{"Rg", 25},
{"Cn", 26},
{"Nh", 27},
{"Fl", 28},
{"Mc", 29},
{"Lv", 30},
{"Ts", 31},
{"Og", 32}
};

const std::map<std::string, int> IsTransMetal
= {
{"H",  0},   
{"He", 0},  
{"Li", 0},  
{"Be", 0},  
{"B",  0},   
{"C",  0},   
{"N",  0},   
{"O",  0},   
{"F",  0},
{"Ne", 0}, 
{"Na", 0}, 
{"Mg", 0}, 
{"Al", 0}, 
{"Si", 0}, 
{"P",  0},  
{"S",  0},  
{"Cl", 0}, 
{"Ar", 0},
{"K",  0},  
{"Ca", 0}, 
{"Sc", 1}, 
{"Ti", 1}, 
{"V",  1},  
{"Cr", 1}, 
{"Mn", 1}, 
{"Fe", 1}, 
{"Co", 1},
{"Ni", 1}, 
{"Cu", 1}, 
{"Zn", 1}, 
{"Ga", 0}, 
{"Ge", 0}, 
{"As", 0}, 
{"Se", 0}, 
{"Br", 0}, 
{"Kr", 0},
{"Rb", 0}, 
{"Sr", 0}, 
{"Y",  1}, 
{"Zr", 1}, 
{"Nb", 1}, 
{"Mo", 1}, 
{"Tc", 1}, 
{"Ru", 1}, 
{"Rh", 1},
{"Pd", 1}, 
{"Ag", 1}, 
{"Cd", 1}, 
{"In", 0}, 
{"Sn", 0}, 
{"Sb", 0}, 
{"Te", 0}, 
{"I",  0},  
{"Xe", 0},
{"Cs", 0}, 
{"Ba", 0}, 
{"La", 1}, 
{"Ce", 1}, 
{"Pr", 1}, 
{"Nd", 1}, 
{"Pm", 1}, 
{"Sm", 1}, 
{"Eu", 1},
{"Gd", 1}, 
{"Tb", 1}, 
{"Dy", 1}, 
{"Ho", 1}, 
{"Er", 1}, 
{"Tm", 1}, 
{"Yb", 1}, 
{"Lu", 1}, 
{"Hf", 1},
{"Ta", 1}, 
{"W",  1},  
{"Re", 1}, 
{"Os", 1}, 
{"Ir", 1}, 
{"Pt", 1}, 
{"Au", 1}, 
{"Hg", 1}, 
{"Tl", 0},
{"Pb", 0}, 
{"Bi", 0}, 
{"Po", 0}, 
{"At", 0}, 
{"Rn", 0}, 
{"Fr", 0},
{"Ra", 0},
{"Ac", 1},
{"Th", 1},
{"Pa", 1},
{"U" , 1},
{"Np", 1},
{"Pu", 1},
{"Am", 1},
{"Cm", 1},
{"Bk", 1},
{"Cf", 1},
{"Es", 1},
{"Fm", 1},
{"Md", 1},
{"No", 1},
{"Lr", 1},
{"Rf", 1},
{"Db", 1},
{"Sg", 1},
{"Bh", 1},
{"Hs", 1},
{"Mt", 1},
{"Ds", 1},
{"Rg", 1},
{"Cn", 1},
{"Nh", 0},
{"Fl", 0},
{"Mc", 0},
{"Lv", 0},
{"Ts", 0},
{"Og", 0}
};

}

#endif

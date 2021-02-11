//
//
//  Copyright (C) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef MINILIB_H
#define MINILIB_H

#include <string>
#include <GraphMol/RDKitBase.h>
#include <memory>
#include "lrucache.hpp"

class JSMol {
 public:
  JSMol() : d_mol(nullptr){};
  JSMol(RDKit::RWMol *mol) : d_mol(mol){};
//  std::string get_smiles() const;
  bool get_substruct_match(const JSMol &q) const;
//  std::string get_morgan_fp(unsigned int radius, unsigned int len) const;
//  std::string get_morgan_fp() const { return get_morgan_fp(2, 2048); };
//  std::string get_sub_fp(unsigned int minPath, unsigned int maxPath, unsigned int fpSize, unsigned int nBitsPerHash) const;
//  std::string get_sub_fp() const { return get_sub_fp(1,6,1024,2); };
  bool is_valid() const { return d_mol.get() != nullptr; };
  std::unique_ptr<RDKit::RWMol> d_mol;
};

JSMol *get_mol(const std::string &input);
JSMol *get_qmol(const std::string &input);
std::string version();
static cache::lru_cache<std::string, JSMol> qmol_cache(3);

#endif
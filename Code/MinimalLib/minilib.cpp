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
#include <string>
#include <iostream>
#include "minilib.h"

#include <RDGeneral/versions.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/Fingerprints/MorganFingerprints.h>
#include <GraphMol/Fingerprints/Fingerprints.h>
#include <DataStructs/BitOps.h>

using namespace RDKit;

namespace {
RWMol *mol_from_input(const std::string &input) {
  RWMol *res = nullptr;
  if (input.find("M  END") != std::string::npos) {
    bool sanitize = false;
    bool removeHs = true;
    bool strictParsing = false;
    res = MolBlockToMol(input, sanitize, removeHs, strictParsing);
  } else {
    SmilesParserParams ps;
    ps.sanitize = false;
    res = SmilesToMol(input, ps);
  }
  if (res) {
    try {
      MolOps::sanitizeMol(*res);
      MolOps::assignStereochemistry(*res, true, true, true);
    } catch (...) {
      delete res;
      res = nullptr;
    }
  }
  return res;
}

RWMol *qmol_from_input(const std::string &input) {
  RWMol *res = nullptr;
  if (input.find("M  END") != std::string::npos) {
    bool sanitize = false;
    bool removeHs = true;
    bool strictParsing = false;
    res = MolBlockToMol(input, sanitize, removeHs, strictParsing);
  } else {
    res = SmartsToMol(input);
  }
  return res;
}
}  // namespace

std::string JSMol::get_smiles() const {
  if (!d_mol) return "";
  return MolToSmiles(*d_mol);
}


bool JSMol::get_substruct_match(const JSMol &q) const {
  if (!d_mol || !q.d_mol) return false;

  MatchVectType match;
  if (SubstructMatch(*d_mol, *(q.d_mol), match)) {
    return true;
  }

  else {
      return false;
  }
}

std::string JSMol::get_sub_fp(unsigned int minPath, unsigned int maxPath, unsigned int fpSize, unsigned int nBitsPerHash) const {
  if (!d_mol) return "";
  auto fp = RDKit::RDKFingerprintMol(*d_mol, minPath, maxPath, fpSize, nBitsPerHash);
  std::string res = BitVectToText(*fp);
  delete fp;
  return res;
}

std::string JSMol::get_morgan_fp(unsigned int radius,
                                 unsigned int fplen) const {
  if (!d_mol) return "";
  auto fp = MorganFingerprints::getFingerprintAsBitVect(*d_mol, radius, fplen);
  std::string res = BitVectToText(*fp);
  delete fp;
  return res;
}

JSMol *get_mol(const std::string &input) {
  RWMol *mol = mol_from_input(input);
  return new JSMol(mol);
}

JSMol *get_qmol(const std::string &input) {
  RWMol *mol = qmol_from_input(input);
  return new JSMol(mol);
}

std::string version() { return std::string(rdkitVersion); }

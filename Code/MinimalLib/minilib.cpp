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
#include <GraphMol/Depictor/RDDepictor.h>
#include <GraphMol/CIPLabeler/CIPLabeler.h>
#include <DataStructs/BitOps.h>

using namespace RDKit;

std::string process_details(const std::string &details, unsigned int &width,
                            unsigned int &height, int &offsetx, int &offsety,
                            std::string &legend, std::vector<int> &atomIds,
                            std::vector<int> &bondIds) {
  rj::Document doc;
  doc.Parse(details.c_str());
  if (!doc.IsObject()) return "Invalid JSON";

  if (doc.HasMember("atoms")) {
    if (!doc["atoms"].IsArray()) {
      return "JSON doesn't contain 'atoms' field, or it is not an array";
    }
    for (const auto &molval : doc["atoms"].GetArray()) {
      if (!molval.IsInt()) return ("Atom IDs should be integers");
      atomIds.push_back(molval.GetInt());
    }
  }
  if (doc.HasMember("bonds")) {
    if (!doc["bonds"].IsArray()) {
      return "JSON contain 'bonds' field, but it is not an array";
    }
    for (const auto &molval : doc["bonds"].GetArray()) {
      if (!molval.IsInt()) return ("Bond IDs should be integers");
      bondIds.push_back(molval.GetInt());
    }
  }

  if (doc.HasMember("width")) {
    if (!doc["width"].IsUint()) {
      return "JSON contains 'width' field, but it is not an unsigned int";
    }
    width = doc["width"].GetUint();
  }

  if (doc.HasMember("height")) {
    if (!doc["height"].IsUint()) {
      return "JSON contains 'height' field, but it is not an unsigned int";
    }
    height = doc["height"].GetUint();
  }

  if (doc.HasMember("offsetx")) {
    if (!doc["offsetx"].IsInt()) {
      return "JSON contains 'offsetx' field, but it is not an int";
    }
    offsetx = doc["offsetx"].GetInt();
  }

  if (doc.HasMember("offsety")) {
    if (!doc["offsety"].IsInt()) {
      return "JSON contains 'offsety' field, but it is not an int";
    }
    offsety = doc["offsety"].GetInt();
  }

  if (doc.HasMember("legend")) {
    if (!doc["legend"].IsString()) {
      return "JSON contains 'legend' field, but it is not a string";
    }
    legend = doc["legend"].GetString();
  }

  return "";
}

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

std::string get_sub_fp(unsigned int minPath, unsigned int maxPath, unsigned int fpSize, unsigned int nBitsPerHash) const {
  if (!d_mol) return "";
  auto fp = Chem.RDKFingerprint(mol, minPath, maxPath, fpSize, nBitsPerHash)
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
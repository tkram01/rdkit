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

const assert = require('assert');
var initRDKitModule = require("../demo/RDKit_minimal.js");
var RDKitModule;

// the goal here isn't to be comprehensive (the RDKit has tests for that),
// just to make sure that the wrappers are working as expected
function test_basics(){
    var bmol = RDKitModule.get_mol("c1ccccc");
    assert.equal(bmol.is_valid(),0);
    
    var mol = RDKitModule.get_mol("c1ccccc1O");
    assert.equal(mol.is_valid(),1);
    assert.equal(mol.get_smiles(),"Oc1ccccc1");

    var fp1 = mol.get_morgan_fp();
    assert.equal(fp1.length,2048);
    assert.equal((fp1.match(/1/g)||[]).length,11);
    var fp2 = mol.get_morgan_fp(0,512);
    assert.equal(fp2.length,512);
    assert.equal((fp2.match(/1/g)||[]).length,3);
    
    var qmol = RDKitModule.get_qmol("Oc(c)c");
    assert.equal(qmol.is_valid(),1);
    var match = mol.get_substruct_match(qmol);
    assert(match, true);
}

initRDKitModule().then(function(instance) {
    RDKitModule = instance;
    console.log(RDKitModule.version());
    test_basics();
    console.log("Tests finished successfully");
});

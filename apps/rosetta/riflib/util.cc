// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <riflib/util.hh>
#include <riflib/scaffold/util.hh>
#include <riflib/rosetta_field.hh>
#include <riflib/scaffold/MultithreadPoseCloner.hh>
#include <numeric/agglomerative_hierarchical_clustering.hh>
#include <scheme/chemical/stub.hh>
#include <scheme/numeric/rand_xform.hh>

#include <ObjexxFCL/format.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <utility/file/file_sys_util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>


#include <core/scoring/rms_util.hh>
#include <core/select/util.hh>

#include <boost/functional/hash/hash.hpp>


namespace devel {
namespace scheme {

void
add_res_if_not_CGP(
	core::pose::Pose const & pose,
	int ir,
	utility::vector1<core::Size> & res,
	bool verbose = false
){
	if( pose.residue(ir).name3() == "CYD" ){ if(verbose) std::cout << "get_res skip CYD" << std::endl; return; }
	if( pose.residue(ir).name3() == "CYS" ){ if(verbose) std::cout << "get_res skip CYS" << std::endl; return; }
	if( pose.residue(ir).name3() == "GLY" ){ if(verbose) std::cout << "get_res skip GLY" << std::endl; return; }
	if( pose.residue(ir).name3() == "PRO" ){ if(verbose) std::cout << "get_res skip PRO" << std::endl; return; }
	res.push_back( ir );
}

utility::vector1<core::Size>
get_res(
	std::string fname,
	core::pose::Pose const & pose,
	bool nocgp
 ){
 	using core::Size;
	utility::vector1<core::Size> res;
	if(fname==""){
		for(Size ir = 1; ir <= pose.size(); ++ir){
			if( nocgp )	add_res_if_not_CGP( pose, ir, res );
			else res.push_back(ir);
		}
		return res;
	}
	runtime_assert_msg( utility::file::file_exists( fname ), "res list file doesn't exist: " + fname );
	utility::io::izstream in(fname);
	std::string s;
	while(in>>s){
		try {
			utility::vector1<std::string> splt = utility::string_split(s,'-');
			// cout << splt << endl;
			if(splt.size()==1){
				// std::cout << "'" << splt.front() << "'" << std::endl;
				Size ir = boost::lexical_cast<Size>(splt.front());
				runtime_assert_msg( ir > 0 && ir <= pose.size(), "residue number out of bounds "+boost::lexical_cast<std::string>(ir) );
				if( nocgp )	add_res_if_not_CGP( pose, ir, res );
				else res.push_back(ir);
			}
			else if(splt.size()==2){
				Size lb = boost::lexical_cast<Size>(splt.front());
				Size ub = boost::lexical_cast<Size>(splt.back());
				if( ub < lb ) std::swap(lb,ub);
				for(Size ir = lb; ir <= ub; ++ir){
					runtime_assert_msg( ir > 0 && ir <= pose.size(), "residue number out of bounds "+boost::lexical_cast<std::string>(ir) );
				if( nocgp )	add_res_if_not_CGP( pose, ir, res );
				else res.push_back(ir);
				}
			} else {
				utility_exit_with_message("get_res: can't parse res line "+s);
			}
		} catch( boost::bad_lexical_cast ){
			utility_exit_with_message("get_res: can't parse res line "+s);
		}
	}
	in.close();

	std::set<core::Size> uniq;
	std::copy(res.begin(), res.end(), std::inserter(uniq, uniq.begin()));
	res.clear();
	std::copy(uniq.begin(), uniq.end(), std::back_inserter(res));
	std::sort( res.begin(), res.end() );
	if( nocgp ){
		for( int iri = 1; iri <= res.size(); ++iri ){
			int ir = res[iri];
			// std::cout << "selected res " << ir << " " << pose.residue(ir).name3() << std::endl;
			runtime_assert( pose.residue(ir).name3() != "CYD" );
			runtime_assert( pose.residue(ir).name3() != "GLY" );
			runtime_assert( pose.residue(ir).name3() != "PRO" );
		}
	}
	// utility_exit_with_message("oriestnd");

	return res;
 }

utility::vector1<core::Size> get_designable_positions_best_guess(
	  core::pose::Pose pose
	, bool noloops
	, bool nocgp
	, bool mutate_to_val
	, bool nohelix
	, bool nostrand
 ){
 	core::scoring::dssp::Dssp dssp( pose );
 	dssp.insert_ss_into_pose( pose );

 	core::pose::Pose allgly = pose;
	core::pose::Pose allval = pose;
 	pose_to_gly(allgly);
 	if (mutate_to_val) pose_to_val(allval);

 	core::scoring::ScoreFunctionOP sf = core::scoring::get_score_function(true);
 	sf->score(allgly);

	core::scoring::methods::EnergyMethodOptions myopt = sf->energy_method_options();
	myopt.hbond_options().decompose_bb_hb_into_pair_energies(true);
	sf->set_energy_method_options(myopt);
 	sf->score(pose);
	core::scoring::hbonds::HBondSet hbset;
	core::scoring::hbonds::fill_hbond_set(pose,false,hbset);

	std::vector<float> sc2bb_energy(pose.size(),false);
	for(core::Size ihb = 1; ihb <= hbset.nhbonds(); ++ihb){
		core::scoring::hbonds::HBond const & hb(hbset.hbond(ihb));
		// noly donot design capping residues???
		if(  hb.don_hatm_is_protein_backbone() && !hb.acc_atm_is_protein_backbone() && pose.secstruct( hb.acc_res() ) == 'L' ) sc2bb_energy[hb.acc_res()-1] += hb.energy();
		if( !hb.don_hatm_is_protein_backbone() &&  hb.acc_atm_is_protein_backbone() && pose.secstruct( hb.don_res() ) == 'L' ) sc2bb_energy[hb.don_res()-1] += hb.energy();
	}

	utility::vector1<core::Size> res;
	// Real calc_per_atom_sasa( pose::Pose const & pose, id::AtomID_Map< Real > & atom_sasa, utility::vector1< Real > & rsd_sasa,
	// Real const probe_radius, bool const use_big_polar_H = false );
	core::id::AtomID_Map< core::Real > atom_sasa;
	utility::vector1< core::Real > rsd_sasa;
	core::scoring::calc_per_atom_sasa( allval, atom_sasa, rsd_sasa, 2.1 );
	for( int ir = 1; ir <= allval.size(); ++ir ){
		// std::cout << pose.secstruct(ir) << std::endl;
		bool isloop = pose.secstruct(ir) == 'L';
		bool ishelix = pose.secstruct(ir) == 'H';
		bool isstrand = pose.secstruct(ir) == 'E';
		int natoms = allval.residue(ir).nheavyatoms()-allval.residue(ir).last_backbone_atom();
		core::Real scsasa = 0;
		for( int ia = allval.residue(ir).first_sidechain_atom(); ia <= allval.residue(ir).natoms(); ++ia ){
			scsasa += atom_sasa[core::id::AtomID(ia,ir)];
		}
		float scsasa_per_atom = scsasa / float(natoms);

		bool is_exposed = scsasa_per_atom > 5.5;
		is_exposed |= scsasa > 18.0;
		// Longxing changed the cut off value below to -0.6, as most of the real capping Hbond score is lower than -1.0.
		bool has_bbsc = sc2bb_energy[ir-1] < -0.6;
		// if(pose.residue(ir).nchi()==4) has_bbsc = false;

		float datr =    pose.energies().residue_total_energies(ir)[core::scoring::fa_atr]
		             -allgly.energies().residue_total_energies(ir)[core::scoring::fa_atr];
		// float dhb  =    pose.energies().residue_total_energies(ir)[core::scoring::hbond_sc]
		//              -allgly.energies().residue_total_energies(ir)[core::scoring::hbond_sc];
		//       dhb +=    pose.energies().residue_total_energies(ir)[core::scoring::fa_elec]
		//              -allgly.energies().residue_total_energies(ir)[core::scoring::fa_elec];

		bool lowenergy = false;
		// lowenergy |= (( datr ) / float(natoms+1)) < -0.9;
		// lowenergy |= ( dhb  ) < -1.8;
		lowenergy |= has_bbsc;

		bool design_res = !lowenergy && is_exposed;

		if(noloops) design_res = design_res & !isloop;
		if(nohelix) design_res = design_res & !ishelix;
		if(nostrand) design_res = design_res & !isstrand;

		if( design_res ){
			if( nocgp )	add_res_if_not_CGP( pose, ir, res );
			else res.push_back(ir);
		} else {
			// std::cout << "pruned res by sasa " << ir << std::endl;
		}
	}
	runtime_assert_msg( res.size() > 0, "no residues selected!" );
	std::cout << "get_designable_positions_best_guess keep " << res.size() << " of " << pose.size() << std::endl;
	std::cout << "select designable_positions_best_guess=resi " << res[1];
	for(int i = 2; i <= res.size(); ++i ) std::cout << "+"<<res[i];
	std::cout << std::endl;

	// utility_exit_with_message("foo");

	return res;
}

std::string
get_res_list_hash( utility::vector1<core::Size> const & reslist ){
	if ( reslist.size() == 0 ) {
		utility_exit_with_message("There were no residues specified on your scaffold!");
	}
	assert( reslist.size() > 0 );
	std::size_t reslist_hash = boost::hash_range( reslist.begin(), reslist.end() );
	std::string hashstr = boost::lexical_cast<std::string>(reslist_hash);
	return hashstr.substr(hashstr.size()-8,8);

}



std::string KMGT(double const & x, int const & w, int const & d){
	using ObjexxFCL::format::F;
	if( x < 1e3  ) return F( w, d, x/1e0  )+" ";
	if( x < 1e6  ) return F( w, d, x/1e3  )+"K";
	if( x < 1e9  ) return F( w, d, x/1e6  )+"M";
	if( x < 1e12 ) return F( w, d, x/1e9  )+"G";
	if( x < 1e15 ) return F( w, d, x/1e12 )+"T";
	if( x < 1e18 ) return F( w, d, x/1e15 )+"P";
	if( x < 1e21 ) return F( w, d, x/1e18 )+"E";
	if( x < 1e24 ) return F( w, d, x/1e21 )+"Z";
	else           return F( w, d, x/1e24 )+"Y";

}

void pose_to_ala( core::pose::Pose & pose ){
	core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
	core::conformation::ResidueOP ala = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map("ALA") );
	for( int ir = 1; ir <= pose.size(); ++ir ){
		if( ! pose.residue(ir).is_protein()   ) continue;
		if(   pose.residue(ir).name3()=="GLY" ) continue;
		if(   pose.residue(ir).name3()=="PRO" ) continue;
		if(   pose.residue(ir).name3()=="CYD" ) continue;
		pose.replace_residue( ir, *ala, true );
	}
}
void pose_to_gly( core::pose::Pose & pose ){
	core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
	core::conformation::ResidueOP gly = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map("GLY") );
	for( int ir = 1; ir <= pose.size(); ++ir ){
		if( ! pose.residue(ir).is_protein()   ) continue;
		if(   pose.residue(ir).name3()=="GLY" ) continue;
		if(   pose.residue(ir).name3()=="PRO" ) continue;
		if(   pose.residue(ir).name3()=="CYD" ) continue;
		pose.replace_residue( ir, *gly, true );
	}
}
void pose_to_ala( core::pose::Pose & pose, utility::vector1<core::Size> const & res_sel ){
	core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
	core::conformation::ResidueOP ala = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map("ALA") );
	for( int iri = 1; iri <= res_sel.size(); ++iri ){
		int ir = res_sel[iri];
		if( ! pose.residue(ir).is_protein()   ) continue;
		if(   pose.residue(ir).name3()=="GLY" ) continue;
		if(   pose.residue(ir).name3()=="PRO" ) continue;
		if(   pose.residue(ir).name3()=="CYD" ) continue;
		pose.replace_residue( ir, *ala, true );
	}
}
void pose_to_val( core::pose::Pose & pose ){
	core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
	core::conformation::ResidueOP val = core::conformation::ResidueFactory::create_residue( rts.lock()->name_map("VAL") );
	for( int ir = 1; ir <= pose.size(); ++ir ){
		if( ! pose.residue(ir).is_protein()   ) continue;
		if(   pose.residue(ir).name3()=="GLY" ) continue;
		if(   pose.residue(ir).name3()=="PRO" ) continue;
		if(   pose.residue(ir).name3()=="CYD" ) continue;
		pose.replace_residue( ir, *val, true );
	}
}

std::string
open_for_read_on_path(
	std::vector<std::string> const & path,
	std::string fname,
	utility::io::izstream & in
){
	for( auto const & dir : path ){
		if( utility::file::file_exists( dir+"/"+fname ) ){
			in.close();
			in.clear();
			in.open( dir+"/"+fname, std::ios::binary );
			if( in.good() ){
				// std::cout << "open_for_read_on_path: opened " << dir + "/" + fname << std::endl;
				return dir+"/"+fname;
			}
		}
	}
	return std::string();
}

std::string
open_for_write_on_path(
	std::vector<std::string> const & path,
	std::string fname,
	utility::io::ozstream & out,
	bool create_directorys /* = false */
){
	for( auto const & dir: path ){
		if( create_directorys && !utility::file::file_exists( dir ) ){
			utility::file::create_directory_recursive( dir );
		}
		if( utility::file::file_exists( dir ) ){
			out.close();
			out.clear();
			out.open( dir+"/"+fname, std::ios::binary );
			if( out.good() ){
				// std::cout << "open_for_write_on_path: opened " << dir + "/" + fname << std::endl;
				return dir+"/"+fname;
			}
		}
	}
	return std::string();

}




void
append_pose_to_pose(
	core::pose::Pose & pose1,
	core::pose::Pose const & pose2,
	bool new_chain
){
	append_subpose_to_pose(pose1, pose2, 1, pose2.size(), new_chain);
}



void
append_subpose_to_pose(
	core::pose::Pose & pose1,
	core::pose::Pose const & pose2,
	core::Size start_res,
	core::Size end_res,
	bool new_chain
){
	if ( pose2.size()<start_res ) {
		std::cerr << "Provided starting residue number " << start_res
			<< " less than number residues in appended pose. Nothing to do." << std::endl;
	}
	pose1.append_residue_by_jump(pose2.residue(start_res), pose1.size() , "", "", new_chain);
	for ( core::Size i=start_res+1; i<=end_res; ++i ) {
		if ( pose2.residue(i).is_lower_terminus() || !pose2.residue(i).is_protein() ) {
			if ( i > 1 && pose2.chain(i) == pose2.chain(i-1) ) {
				pose1.append_residue_by_jump(pose2.residue(i), pose1.size(), "","", false);
			} else {
				pose1.append_residue_by_jump(pose2.residue(i), pose1.size(), "","", true);
			}
		} else {
			pose1.append_residue_by_bond(pose2.residue(i));
		}
	}
}




std::vector<int>
get_rif_atype_map()
{
	core::chemical::AtomTypeSetCOP ats = core::chemical::ChemicalManager::get_instance()->atom_type_set("fa_standard");
	// 147 CH0
	// 148 HS
	// 149 NtrR
	// 150 SH1


	std::vector<int> atypemap;
	atypemap.resize( ats->n_atomtypes()+1, 12345 );
	atypemap[ ats->atom_type_index("CNH2") ] =  1;
	atypemap[ ats->atom_type_index("COO" ) ] =  2;
	atypemap[ ats->atom_type_index("CH1" ) ] =  3;
	atypemap[ ats->atom_type_index("CH2" ) ] =  4;
	atypemap[ ats->atom_type_index("CH3" ) ] =  5;
	atypemap[ ats->atom_type_index("aroC") ] =  6;
	atypemap[ ats->atom_type_index("Ntrp") ] =  7;
	atypemap[ ats->atom_type_index("Nhis") ] =  8;
	atypemap[ ats->atom_type_index("NH2O") ] =  9;
	atypemap[ ats->atom_type_index("Nlys") ] = 10;
	atypemap[ ats->atom_type_index("Narg") ] = 11;
	// atypemap[ ats->atom_type_index("Narg") ] = 6; // hack hack hack hack
	atypemap[ ats->atom_type_index("Npro") ] = 12;
	atypemap[ ats->atom_type_index("OH"  ) ] = 13;
	atypemap[ ats->atom_type_index("ONH2") ] = 14;
	atypemap[ ats->atom_type_index("OOC" ) ] = 15;
	atypemap[ ats->atom_type_index("Oaro") ] = 16;
	atypemap[ ats->atom_type_index("S"   ) ] = 17;
	atypemap[ ats->atom_type_index("Nbb" ) ] = 18;
	atypemap[ ats->atom_type_index("CAbb") ] = 19;
	atypemap[ ats->atom_type_index("CObb") ] = 20;
	atypemap[ ats->atom_type_index("OCbb") ] = 21;

	// -beta stuff
	if(ats->has_atom("CH0" )) atypemap[ ats->atom_type_index("CH0" ) ] = 6;
	// if(ats->has_atom("HS")) atypemap[ ats->atom_type_index("HS"  ) ] = ??;
	if(ats->has_atom("NtrR")) atypemap[ ats->atom_type_index("NtrR") ] = 11;
	if(ats->has_atom("SH1" )) atypemap[ ats->atom_type_index("SH1" ) ] = 17;

        // These are in DNA. Blame Robert Pecoraro if you don't like the typing -bcov
	if(ats->has_atom("Oet2" )) atypemap[ ats->atom_type_index("Oet2" ) ] = 21; //OCbb

	return atypemap;
}

// I was going to use the atom_properties.txt, but turns out it was wrong
std::vector<bool> get_rif_hydrophobic_map() {
	std::vector<int> atypemap = get_rif_atype_map();

	std::vector<bool> hydrophobic_map( atypemap.size(), false );

	hydrophobic_map[ 1] = false; //"CNH2"
	hydrophobic_map[ 2] = false; //"COO" 
	hydrophobic_map[ 3] =  true; //"CH1" 
	hydrophobic_map[ 4] =  true; //"CH2" 
	hydrophobic_map[ 5] =  true; //"CH3" 
	hydrophobic_map[ 6] =  true; //"aroC"
	hydrophobic_map[ 7] = false; //"Ntrp"
	hydrophobic_map[ 8] = false; //"Nhis"
	hydrophobic_map[ 9] = false; //"NH2O"
	hydrophobic_map[10] = false; //"Nlys"
	hydrophobic_map[11] = false; //"Narg"
	hydrophobic_map[12] = false; //"Npro"
	hydrophobic_map[13] = false; //"OH"  
	hydrophobic_map[14] = false; //"ONH2"
	hydrophobic_map[15] = false; //"OOC" 
	hydrophobic_map[16] = false; //"Oaro"
	hydrophobic_map[17] =  true; //"S"   
	hydrophobic_map[18] = false; //"Nbb" 
	hydrophobic_map[19] =  true; //"CAbb"
	hydrophobic_map[20] = false; //"CObb"
	hydrophobic_map[21] = false; //"OCbb"

	return hydrophobic_map;
}

std::vector<bool> get_rif_pi_map() {
	std::vector<int> atypemap = get_rif_atype_map();

	std::vector<bool> hydrophobic_map( atypemap.size(), false );

	hydrophobic_map[ 1] = false; //"CNH2"
	hydrophobic_map[ 2] =  true; //"COO" 
	hydrophobic_map[ 3] = false; //"CH1" 
	hydrophobic_map[ 4] = false; //"CH2" 
	hydrophobic_map[ 5] = false; //"CH3" 
	hydrophobic_map[ 6] =  true; //"aroC"
	hydrophobic_map[ 7] =  true; //"Ntrp"
	hydrophobic_map[ 8] =  true; //"Nhis"
	hydrophobic_map[ 9] = false; //"NH2O"
	hydrophobic_map[10] = false; //"Nlys"
	hydrophobic_map[11] = false; //"Narg"
	hydrophobic_map[12] = false; //"Npro"
	hydrophobic_map[13] = false; //"OH"  
	hydrophobic_map[14] = false; //"ONH2"
	hydrophobic_map[15] = false; //"OOC" 
	hydrophobic_map[16] = false; //"Oaro"
	hydrophobic_map[17] = false; //"S"   
	hydrophobic_map[18] = false; //"Nbb" 
	hydrophobic_map[19] = false; //"CAbb"
	hydrophobic_map[20] = false; //"CObb"
	hydrophobic_map[21] = false; //"OCbb"

	return hydrophobic_map;
}



std::map< core::id::AtomID, core::id::AtomID >
residue_subset_to_CA_atom_map( 
    core::select::residue_selector::ResidueSubset const & subset, 
    core::pose::Pose const & pose ) {

    std::map< core::id::AtomID, core::id::AtomID > atom_map;
    utility::vector1<core::Size> seq_poss = core::select::get_residues_from_subset( subset );
    for ( core::Size seq_pos : seq_poss ) {
        core::conformation::Residue const & res( pose.residue( seq_pos ) );
        if ( res.has( "CA" ) ) {
            core::id::AtomID const id( res.atom_index( "CA" ), seq_pos );
            atom_map[ id ] = id;
        }
    }
    return atom_map;
}



core::Real
subset_CA_rmsd(
    core::pose::Pose const & pose1,
    core::pose::Pose const & pose2,
    core::select::residue_selector::ResidueSubset const & subset,
    bool superimpose) {

    utility::vector1<core::Size> calc_rms_res = core::select::get_residues_from_subset( subset );
    std::map< core::id::AtomID, core::id::AtomID > atom_map = 
        residue_subset_to_CA_atom_map( subset, pose1 );

    if ( superimpose ) {
        return core::scoring::rms_at_corresponding_atoms( pose1, pose2, atom_map, calc_rms_res );
    } else {
        return core::scoring::rms_at_corresponding_atoms_no_super( pose1, pose2, atom_map, calc_rms_res );
    }
}


/// Your poses must be identical. Very important
EigenXform
find_xform_from_identical_pose_to_pose( 
	core::pose::Pose const & pose1,
	core::pose::Pose const & pose2,
	float align_error /* = 0.2 */ ) {

	runtime_assert( pose1.size() == pose2.size() );

	// we can use anything here because they're the same
	utility::vector1<core::Size> scaffold_res;
	get_default_scaffold_res( pose1, scaffold_res );

	core::pose::Pose to_move = pose1;
    ::devel::scheme::pose_to_ala( to_move );
    Eigen::Vector3f to_move_center = pose_center(to_move, scaffold_res);

    core::pose::Pose match_this = pose2;
    ::devel::scheme::pose_to_ala( match_this );
    Eigen::Vector3f match_center = pose_center(match_this,scaffold_res);


   	// prepare all test point now

    // Test original residue 1
    utility::vector1<core::Size> target_res {1};
    std::vector< ::scheme::actor::Atom< Eigen::Vector3f > > to_move_atoms;
    std::vector< ::scheme::actor::Atom< Eigen::Vector3f > > match_atoms;
    
    devel::scheme::get_scheme_atoms( to_move, target_res, to_move_atoms, true );
    devel::scheme::get_scheme_atoms( match_this, target_res, match_atoms, true );


    // Test original midpoint
    core::Size test_res = to_move.size() / 2;

    utility::vector1<core::Size> test_target_res {test_res};
    std::vector< ::scheme::actor::Atom< Eigen::Vector3f > > test_to_move_atoms;
    std::vector< ::scheme::actor::Atom< Eigen::Vector3f > > test_match_atoms;

    devel::scheme::get_scheme_atoms( to_move, test_target_res, test_to_move_atoms, true );
    devel::scheme::get_scheme_atoms( match_this, test_target_res, test_match_atoms, true );


    // Test centered residue 1

    EigenXform to_move_2_center = EigenXform::Identity();
    to_move_2_center.translation() = -to_move_center;
    apply_xform_to_pose( to_move, to_move_2_center );

    std::vector< ::scheme::actor::Atom< Eigen::Vector3f > > to_move_centered_atoms;
    devel::scheme::get_scheme_atoms( to_move, target_res, to_move_centered_atoms, true );


    // done making test points




    EigenXform match_x = ::scheme::chemical::make_stub<EigenXform>(
                                                                match_atoms[0].position(),
                                                                match_atoms[1].position(),
                                                                match_atoms[2].position());
    EigenXform to_move_x = ::scheme::chemical::make_stub<EigenXform>(
                                                                to_move_centered_atoms[0].position(),
                                                                to_move_centered_atoms[1].position(),
                                                                to_move_centered_atoms[2].position());

    EigenXform to_move_centered_2_match = match_x * to_move_x.inverse();
    to_move_centered_2_match.translation() = match_center;

    double centered_error = (to_move_centered_2_match * to_move_centered_atoms[0].position() - match_atoms[0].position()).norm();
	std::cout << "Centered res1 Alignment error :" << centered_error << std::endl;
    runtime_assert( centered_error < align_error );

    EigenXform to_move_2_match = to_move_centered_2_match * to_move_2_center;

    double error = (to_move_2_match * to_move_atoms[0].position() - match_atoms[0].position()).norm();
	std::cout << "res1 Alignment error :" << error << std::endl;
    runtime_assert( error < align_error );


    double test_error = (to_move_2_match * test_to_move_atoms[0].position() - test_match_atoms[0].position()).norm();
    std::cout << "Test res" << test_res << " Alignment error :" << error << std::endl;
    runtime_assert( test_error < align_error );

    return to_move_2_match;
	
}



void
apply_xform_to_pose( core::pose::Pose & pose, EigenXform const & xform) {
	numeric::xyzTransform<float> transform = eigen2xyz( xform );
	pose.apply_transform_Rx_plus_v( transform.R, transform.t );

}

void
apply_xform_to_residue( core::conformation::Residue & residue, EigenXform const & xform) {
	numeric::xyzTransform<float> transform = eigen2xyz( xform );
	residue.apply_transform_Rx_plus_v( transform.R, transform.t );

}

Eigen::Vector3f
pose_center(
    core::pose::Pose const & pose,
    utility::vector1<core::Size> const & useres /* = utility::vector1<core::Size>()*/
){
    typedef numeric::xyzVector<core::Real> Vec;
    Vec cen(0,0,0);
    int count = 0;
    for( int ir = 1; ir <= pose.size(); ++ir ) {
        if( useres.size()==0 || std::find(useres.begin(),useres.end(),ir)!=useres.end() ){
            for( int ia = 1; ia <= pose.residue(ir).nheavyatoms(); ++ia ){
                cen += pose.xyz(core::id::AtomID(ia,ir));
                ++count;
            }
        // } else {
            // std::cout << "pose_center skip " << ir << std::endl;
        }
    }
    cen /= double(count);
    // ::scheme::util::SimpleArray<3,float> center;
    Eigen::Vector3f center;
    center[0] = cen[0];
    center[1] = cen[1];
    center[2] = cen[2];
    return center;
}

void
get_rg_radius(
    core::pose::Pose const & pose,
    float & rg,
    float & radius,
    utility::vector1<core::Size> const & useres /*= utility::vector1<core::Size>()*/,
    bool allatom /*= false*/
){
    Eigen::Vector3f centmp = pose_center( pose, useres );
    numeric::xyzVector<double> cen;
    float maxdis = -9e9, avgdis2 = 0.0;
    for( int i = 0; i < 3; ++i ) cen[i] = centmp[i];
    for( int iri = 1; iri <= useres.size(); ++iri ){
        int ir = useres[iri];
        if( allatom ){
            for( int ia = 1; ia <= pose.residue(ir).nheavyatoms(); ++ia ){
                numeric::xyzVector<double> coord = pose.residue(ir).xyz(ia);
                avgdis2 += cen.distance_squared( coord );
                maxdis = std::max( maxdis, (float)cen.distance( coord ) );
            }
        } else {
            numeric::xyzVector<double> coord;
            if(      pose.residue(ir).has("CB") ) coord = pose.residue(ir).xyz("CB");
            else if( pose.residue(ir).has("CA") ) coord = pose.residue(ir).xyz("CA");
            else                                  coord = pose.residue(ir).nbr_atom_xyz();
            avgdis2 += cen.distance_squared( coord );
            maxdis = std::max( maxdis, (float)cen.distance( coord ) );
        }
    }
    avgdis2 /= useres.size();
    rg = sqrt( avgdis2 );
    radius = maxdis;
}

// in radians
float
angle_between_vectors(
	Eigen::Vector3f const & vec1,
	Eigen::Vector3f const & vec2
) {
	const float dot = vec1.dot(vec2);
	const float norm_product = vec1.norm() * vec2.norm();
	return std::acos( dot / norm_product );
}

void 
xform_pose( 
	core::pose::Pose & pose, 
	numeric::xyzTransform<float> s, 
	core::Size sres/*=1*/, 
	core::Size eres/*=0*/
) {
  if(eres==0) eres = pose.size();
  for(core::Size ir = sres; ir <= eres; ++ir) {
    for(core::Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
      core::id::AtomID const aid(core::id::AtomID(ia,ir));
      pose.set_xyz( aid, s*pose.xyz(aid) );
    }
  }
}



}
}



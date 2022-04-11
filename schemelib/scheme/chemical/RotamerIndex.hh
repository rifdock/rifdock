#ifndef INCLUDED_chemical_RotamerIndex_HH
#define INCLUDED_chemical_RotamerIndex_HH

#include "scheme/util/assert.hh"
#include "scheme/util/str.hh"
#include "scheme/chemical/AtomData.hh"
#include "scheme/chemical/stub.hh"
#include "scheme/chemical/HBondRay.hh"
#include "scheme/io/dump_pdb_atom.hh"
#include "scheme/actor/BackboneActor.hh"

#include <boost/functional/hash.hpp>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzTransform.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/lkball/LK_BallInfo.hh>
#include <core/pose/Pose.hh>

#include <cstdlib>
#include <string>
#include <iostream>
#include <map>
#include <cmath>
#include <Eigen/Dense>

namespace scheme { namespace chemical {


template< class _AtomData >
struct ChemicalIndex {
	typedef _AtomData AtomData;
	std::vector<std::string> resnames_;
	std::map<std::string,int> resname2num_;
	std::vector< std::vector< AtomData > > atomdata_;
	AtomData null_atomdata_;

	ChemicalIndex() : null_atomdata_() {
	}

	bool have_res( std::string resn ) const {
		return resname2num_.find(resn) != resname2num_.end();
	}

	void add_res( std::string resn ){
		if( have_res(resn) ) return;
		resname2num_[resn] = resnames_.size();
		resnames_.push_back(resn);
	}

	void add_atomdata( int restype, int atomnum, AtomData const & data ){
		if( atomdata_.size() <= restype ) atomdata_.resize(restype+1);
		if( atomdata_[restype].size() <= atomnum ) atomdata_[restype].resize(atomnum+1);
		atomdata_[restype][atomnum] = data;
	}
	void add_atomdata( std::string resname, int atomnum, AtomData const & data ){
		add_res( resname );
		add_atomdata( resname2num_[resname], atomnum, data );
	}

	AtomData       & atom_data( int restype, int atomnum )       { if(restype<0) return null_atomdata_; return atomdata_.at(restype).at(atomnum); }
	AtomData const & atom_data( int restype, int atomnum ) const { if(restype<0) return null_atomdata_; return atomdata_.at(restype).at(atomnum); }
	AtomData       & atom_data( std::string const & resname, int atomnum )       { return atomdata_.at( resname2num_.find(resname)->second ).at(atomnum); }
	AtomData const & atom_data( std::string const & resname, int atomnum ) const { return atomdata_.at( resname2num_.find(resname)->second ).at(atomnum); }

	int resname2num( std::string const & resn ) const {
		std::map<std::string,int>::const_iterator i = resname2num_.find(resn);
		if( i == resname2num_.end() ) return -1;
		return i->second;
	}

	void clear(){
		resnames_.clear();
		resname2num_.clear();
		atomdata_.clear();
	}

	bool operator==(ChemicalIndex<AtomData> const & other) const {
		return (
				resnames_      == other.resnames_    &&
				resname2num_   == other.resname2num_ &&
				atomdata_      == other.atomdata_    &&
				null_atomdata_ == other.null_atomdata_
			);
	}

};

namespace impl {

template< class _Atom >
struct Rotamer {
	typedef _Atom Atom;
	// typedef struct { typename Atom::Position position; int32_t type; } AtomPT;
	std::string resname_;
	size_t n_proton_chi_;
	std::vector<float> chi_;
	std::vector<Atom> atoms_;
	std::vector<std::pair<Atom,Atom> > hbonders_;
	std::vector<HBondRay> donors_;
	std::vector<HBondRay> acceptors_;
	int nheavyatoms;
	uint64_t validation_hash() const {
		uint64_t h = (uint64_t)(boost::hash<std::string>()(resname_));
		h ^= (uint64_t)boost::hash<size_t>()(n_proton_chi_);
		for(int i = 0; i < chi_.size(); ++i){
			h ^= (uint64_t)boost::hash<float>()(chi_.at(i));
		}
		return h;
	}
	bool operator==(Rotamer<Atom> const & other) const {
		bool atoms_close = atoms_.size() == other.atoms_.size();
		if( atoms_close ){
			for( int i = 0; i < atoms_.size(); ++i ){
				float dist = (atoms_.at(i).position() - other.atoms_.at(i).position()).norm();
				atoms_close &= dist < 0.00001;
				atoms_close &= atoms_.at(i).type() == other.atoms_.at(i).type();
			}
		}
		return (
				resname_ == other.resname_ &&
				n_proton_chi_ == other.n_proton_chi_ &&
				chi_ == other.chi_ &&
				atoms_close &&
				// atoms_ == other.atoms_ &&
				// hbonders_ == other.hbonders_ &&
				donors_ == other.donors_ &&
				acceptors_ == other.acceptors_ &&
				nheavyatoms == other.nheavyatoms
			);
	}
};

	template<class F=float>
	bool angle_is_close( F a, F b, F d ){
		F diff = std::min( std::abs(a-b+360.0f), std::abs(a-b-360.0f) );
		diff = std::min( diff, std::abs(a-b) );
		return diff <= d;
	}

}

//Adding Specific Rotamers
struct
RotamerSpec{
	std::string resname_;
	std::vector<float> chi_;
	int parent_key_, n_proton_chi_;
};




inline
void
residue_to_identity( core::conformation::ResidueOP & resop ) {
	typedef ::Eigen::Transform<float,3,Eigen::AffineCompact> EigenXform;

	Eigen::Vector3f N ( resop->xyz("N" ).x(), resop->xyz("N" ).y(), resop->xyz("N" ).z() );
	Eigen::Vector3f CA( resop->xyz("CA").x(), resop->xyz("CA").y(), resop->xyz("CA").z() );
	Eigen::Vector3f C ( resop->xyz("C" ).x(), resop->xyz("C" ).y(), resop->xyz("C" ).z() );
	// typedef Eigen::Transform<float,3,Eigen::AffineCompact> EigenXform;
	::scheme::actor::BackboneActor<EigenXform> bbactor( N, CA , C );

	EigenXform xform = bbactor.position().inverse();

	::numeric::xyzMatrix<float> m;
	for(int i = 0; i < 3; ++i){
		for(int j = 0; j < 3; ++j){
			m(i+1,j+1) = xform.rotation()(i,j);
		}
	}
	::numeric::xyzTransform<float> x(m);
	x.t[0] = xform.translation()[0];
	x.t[1] = xform.translation()[1];
	x.t[2] = xform.translation()[2];

	resop->apply_transform_Rx_plus_v( x.R, x.t );
}

inline
core::conformation::ResidueOP
get_residue_at_identity( 
	core::chemical::ResidueType const & rtype,
	std::vector<float> const & chi 
) {

	core::conformation::ResidueOP resop = core::conformation::ResidueFactory::create_residue( rtype );
	runtime_assert( chi.size() == resop->nchi() );
	for(int i = 0; i < chi.size(); ++i){
		resop->set_chi( i+1, chi[i] );
	}

	residue_to_identity( resop );

	return resop;

}

inline
void get_residue_rotspec_params(
	core::conformation::Residue const & res,
	std::string & resn,
	std::vector<float> & chis,
	int & n_proton_chi,
	int & parent_key
) {
	resn = res.name3();
	parent_key = -1;
	n_proton_chi = 0;
	if (resn == "CYS" || resn == "SER" || resn == "THR" || resn == "TYR"){
		n_proton_chi = 1;
	}
	chis.clear();
	chis.resize(res.nchi());
	for (int n_chi = 0; n_chi < chis.size(); n_chi++){
		chis.at(n_chi) = res.chi(n_chi+1);
	}
}



//template<class _Atom, class RotamerGenerator, class Xform>
struct
RotamerIndexSpec
{
	//typedef _Atom Atom;
	//typedef impl::Rotamer<Atom> Rotamer;
	// typedef RotamerIndexSpec THIS
	// RotamerIndexSpec(){
	// 	this -> std::vector<RotamerSpec> rot_specs;
	// }
	std::vector<RotamerSpec> rot_specs;

	
	size_t size() const { return rot_specs.size(); }
	void clear() {rot_specs.clear();}
	std::string resname(size_t i) const { return rot_specs.at(i).resname_; }
	scheme::chemical::RotamerSpec get_rotspec(int i) const {return rot_specs.at(i);}


	int get_matching_rot(
		std::string & resn,
		std::vector<float> & chis,
		int & n_proton_chi,
		float tolerance=5.0f
	) {
		for( int irot = 0; irot < size(); ++irot ){

			if( resname(irot) != resn ) continue;
			bool match = true;
			for (int n_chi = 0; n_chi < get_rotspec(irot).chi_.size(); ++n_chi){
				match &= ::scheme::chemical::impl::angle_is_close( get_rotspec(irot).chi_.at(n_chi), chis.at(n_chi), tolerance );
			}

			if (match)
				return irot;				
		}// end loop over all existing rot_spec rotamers
		return -1;
	}

	// returns irot
	int get_matching_rot( core::conformation::Residue const & res, float tolerance=5.0f) {
		std::string resn;
		std::vector<float> chis;
		int n_proton_chi;
		int parent_key;
		get_residue_rotspec_params( res, resn, chis, n_proton_chi, parent_key );

		return get_matching_rot( resn, chis, n_proton_chi, tolerance );
	}

	int add_rotamer(std::string resname, std::vector<float> const & chi, int n_proton_chi, int parent_key = -1){
		RotamerSpec S;
		S.resname_ = resname;
		S.chi_ = chi;
		S.n_proton_chi_ = n_proton_chi;
		S.parent_key_ = parent_key;
		rot_specs.push_back(S);
		return rot_specs.size()-1;
	};
	template<class T>
	void fill_rotamer_index(T& rot_index) const {
		//std::cout << rot_index.size() << std::endl;
		//std::cout << rot_specs.size() << std::endl;
		for(auto rot_spec : rot_specs){
			//std::cout << rot_spec.resname_ << std::endl;
			if (rot_index.d_l_map_.find(rot_spec.resname_) == rot_index.d_l_map_.end()) {
				rot_index.add_rotamer(rot_spec.resname_,rot_spec.chi_,rot_spec.n_proton_chi_,rot_spec.parent_key_);
			} else {
				rot_index.add_rotamer(rot_spec.resname_,rot_spec.chi_,rot_spec.n_proton_chi_,rot_spec.parent_key_, true);
			}
		}
		std::cout << "start building rotamer index from fill rotamers" << std::endl;
		rot_index.build_index();
		std::cout << "finish build_index..." << std::endl;
	}
	void save(std::ostream &out){
		//out << "test save method";
		for (int i = 0; i < rot_specs.size();i++){
			out << rot_specs.at(i).resname_ << "\t";
			for (int k = 0; k < rot_specs.at(i).chi_.size(); k++){
				out << rot_specs.at(i).chi_.at(k) << "\t";	
			}
			out << rot_specs.at(i).n_proton_chi_ << "\t";
			out << rot_specs.at(i).parent_key_ << "\t";
			out << "\n";
		}

	}
	void load(std::istream &in){
		std::string rot_line;
		while (std::getline(in,rot_line)){
			std::istringstream iss(rot_line);
			std::vector<std::string> rot_info((std::istream_iterator<std::string>(iss)),std::istream_iterator<std::string>());
			
			//store the rot_info
			std::string my_resn = rot_info.at(0);
			int my_parent_key = std::stoi(rot_info.at(rot_info.size()-1),nullptr);
			int my_n_proton_chi = std::stoi(rot_info.at(rot_info.size()-2),nullptr);
			std::vector<float> my_chi;

			for (std::vector<std::string>::iterator it = rot_info.begin()+1;it != rot_info.end()-2;++it){
				my_chi.push_back(strtof((*it).c_str(),nullptr));

			}
			add_rotamer(my_resn,my_chi,my_n_proton_chi,my_parent_key);
			//quick test;
			// for (int i = 0; i < rot_info.size(); i++){
			//  	std::cout << rot_info[i] << "\t";
			// }
			// std::cout << "\n";
			// std::cout << "my_prot_chi: " << my_n_proton_chi << std::endl;
			// std::cout << "my_parent_key: " << my_parent_key << std::endl;
			// for (int i = 0; i < my_chi.size(); i++){
			// 	std::cout << my_chi[i] << "\t" << std::endl;
			// }
		}

		//rot_specs
		//clear();
		// size_t n;
		// in.read((char*)&n,sizeof(size_t));
		// // std::cout << "read: " << n << std::endl;
		// for( size_t i = 0; i < n; ++i){
		// 	uint8_t nresn;
		// 	in.read((char*)&nresn,sizeof(uint8_t));
		// 	std::string resn(nresn,0);
		// 	for( int k = 0; k < nresn; ++k ){
		// 		in.read((char*)&resn[k],sizeof(char));
		// 	}
		// 	// std::cout << "read " << i << " " << resn << std::endl;
		// 	uint8_t n_chi;
		// 	in.read((char*)&n_chi,sizeof(uint8_t));
		// 	std::vector<float> mychi(n_chi);
		// 	for( int k = 0; k < n_chi; ++k ){
		// 		in.read((char*)&mychi[k],sizeof(float));
		// 	}
		// 	// std::cout << "read " << i << " " << mychi.size() << " " << n_chi << std::endl;
		// 	int npchi;
		// 	in.read((char*)&npchi,sizeof(int));
		// 	// std::cout << "read " << i << " " << npchi << std::endl;
		// 	int parent;
		// 	in.read((char*)&parent,sizeof(int));
		// 	// std::cout << "read " << i << " " << parent << std::endl;
		// 	add_rotamer(resn, mychi, npchi, parent);
		// }
		// std::string endtag;
		// in >> endtag;
		// ALWAYS_ASSERT(endtag=="load_rot_index_end");
		//rot_specs.build_index();
	}

	// This is not a fast function!! Used for dumping the rotamer spec
	core::conformation::ResidueOP
	get_rotamer_at_identity( size_t irot, std::string lresname = "" ) const {
		core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
		core::conformation::ResidueOP rop;
		if (lresname == "") {
			rop = get_residue_at_identity( rts.lock()->name_map( this->resname(irot)), get_rotspec(irot).chi_ );
		} else {
			core::chemical::ResidueTypeCOP rt = rts.lock() -> name_mapOP( lresname );
			core::chemical::ResidueTypeCOP d_rt = rts.lock() -> get_d_equivalent(rt);
			rop = get_residue_at_identity( *d_rt, get_rotspec(irot).chi_ );
		}
		return rop;
	}

	void
	dump_all_rotspec_rotamers( std::string const & filename ) const {
		::numeric::xyzVector< core::Real > x_unit( 1, 0, 0 );
		::numeric::xyzMatrix< core::Real > identity = ::numeric::xyzMatrix< core::Real >::identity();

		core::Real const X_SCALE = 5;

		core::pose::Pose pose;
		for ( size_t irot = 0; irot < this->size(); irot++ ) {
			core::conformation::ResidueOP res = get_rotamer_at_identity( irot );
			res->apply_transform_Rx_plus_v( identity, x_unit * irot * X_SCALE );
			pose.append_residue_by_jump( *res, 1 );
		}

		pose.dump_pdb( filename );
	}

};

template<class _Atom, class RotamerGenerator, class Xform>
struct RotamerIndex {
	BOOST_STATIC_ASSERT(( ( boost::is_same< _Atom, typename RotamerGenerator::Atom >::value ) ));
	typedef RotamerIndex<_Atom,RotamerGenerator,Xform> THIS;
	typedef _Atom Atom;
	typedef impl::Rotamer<Atom> Rotamer;

	size_t size() const { return rotamers_.size(); }
	size_t n_primary_rotamers() const { return n_primary_rotamers_; }
	std::string resname(size_t i) const { return rotamers_.at(i).resname_; }
	std::string oneletter(size_t i) const { 
		if ( oneletter_map_.find(rotamers_.at(i).resname_) == oneletter_map_.end() ) {
			return rotamers_.at(i).resname_;
		} else {
			return oneletter_map_.at(rotamers_.at(i).resname_);
		}
	}
	size_t natoms( size_t i ) const { return rotamers_.at(i).atoms_.size(); }
	size_t nheavyatoms( size_t i ) const { return rotamers_.at(i).nheavyatoms; }
	size_t nchi( size_t i ) const { return rotamers_.at(i).chi_.size(); }
	size_t nchi_noproton( size_t i ) const { return nchi(i)-nprotonchi(i); }
	size_t nprotonchi( size_t i ) const { return rotamers_.at(i).n_proton_chi_; }
	float   chi( size_t i, size_t j ) const { return rotamers_.at(i).chi_.at(j); }
	size_t nhbonds( size_t i ) const { return rotamers_.at(i).hbonders_.size(); }
	Atom const & atom( size_t i, size_t iatm) const { return rotamers_.at(i).atoms_.at(iatm); }
	Atom const & hbond_atom1( size_t i, size_t ihb ) const { return rotamers_.at(i).hbonders_.at(ihb).first; }
	Atom const & hbond_atom2( size_t i, size_t ihb ) const { return rotamers_.at(i).hbonders_.at(ihb).second; }
	size_t parent_irot( size_t i ) const { return parent_rotamer_.at(i); }
	int same_struct_start_chi(size_t irot) const { return (resname(irot)=="ILE" || resname(irot) == "DIL") ? 1 : 2; }
	Rotamer const & rotamer( int irot ) const { return rotamers_.at(irot); }

	std::map<std::string,std::string> oneletter_map_;
	// map between d and l version of aa
	std::map<std::string,std::string> d_l_map_;
	std::vector<bool> is_d_;

	ChemicalIndex<AtomData> chem_index_;

	std::vector< Rotamer > rotamers_;

	int n_primary_rotamers_ = 0;
	bool seen_child_rotamer_ = false;
	std::vector< int >t_rotamer_;
	RotamerGenerator rotgen_;

	std::vector<int> parent_rotamer_;

	typedef std::map<std::string,std::pair<int,int> > BoundsMap;
	BoundsMap bounds_map_;

	typedef std::vector< std::pair<int,int> > ChildMap; // index is "primary" index, not rotamer number
	ChildMap child_map_;
	std::vector<int> protonchi_parent_of_;
	int ala_rot_ = -1;

	std::vector<int> structural_parents_;
	std::vector<int> structural_parent_of_;
	std::vector<Xform> to_structural_parent_frame_;

	// this is egregious and is only used with GridScoring
	// per_thread_rotamers_[thread][irot];
	// TODO: probably broke it with d_aa code need to be fix
	std::vector<std::vector<core::conformation::ResidueOP>> per_thread_rotamers_;
	std::vector<std::vector<core::scoring::lkball::LKB_ResidueInfoOP>> per_thread_lkbrinfo_;


	RotamerIndex(){
		this->fill_oneletter_map( oneletter_map_ );
		this->fill_d_l_map( d_l_map_ );
	}

	void clear(){
		n_primary_rotamers_ = 0;
		seen_child_rotamer_ = false;
		ala_rot_ = -1;
		chem_index_.clear();
		rotamers_.clear();
		parent_rotamer_.clear();
		bounds_map_.clear();
		child_map_.clear();
		protonchi_parent_of_.clear();
		structural_parents_.clear();
		structural_parent_of_.clear();
		to_structural_parent_frame_.clear();
		// keep oneletter_map_
	}

	void save(std::ostream & ostrm) const {
		size_t n = size();
		ostrm.write((char*)&n,sizeof(size_t));
		for( size_t i = 0; i < n; ++i ){
			uint8_t nresn = resname(i).size();
			ostrm.write((char*)&nresn,sizeof(uint8_t));
			for( int k = 0; k < nresn; ++k ){
				ostrm.write((char*)&resname(i)[k],sizeof(char));
			}
			uint8_t n_chi = nchi(i);
			ostrm.write((char*)&n_chi,sizeof(uint8_t));
			for( int k = 0; k < n_chi; ++k ){
				float tmpchi = chi(i,k);
				ostrm.write((char*)&tmpchi,sizeof(float));
			}
			int npchi = nprotonchi(i);
			ostrm.write((char*)&npchi,sizeof(int));
			int parent = parent_rotamer_.at(i);
			if( parent == i ) parent = -1;
			ostrm.write((char*)&parent,sizeof(int));
		}
		ostrm << "rot_index_end";
	}

	void load(std::istream & istrm){
		clear();
		size_t n;
		istrm.read((char*)&n,sizeof(size_t));
		// std::cout << "read: " << n << std::endl;
		for( size_t i = 0; i < n; ++i){
			uint8_t nresn;
			istrm.read((char*)&nresn,sizeof(uint8_t));
			std::string resn(nresn,0);
			for( int k = 0; k < nresn; ++k ){
				istrm.read((char*)&resn[k],sizeof(char));
			}
			// std::cout << "read " << i << " " << resn << std::endl;
			uint8_t n_chi;
			istrm.read((char*)&n_chi,sizeof(uint8_t));
			std::vector<float> mychi(n_chi);
			for( int k = 0; k < n_chi; ++k ){
				istrm.read((char*)&mychi[k],sizeof(float));
			}
			// std::cout << "read " << i << " " << mychi.size() << " " << n_chi << std::endl;
			int npchi;
			istrm.read((char*)&npchi,sizeof(int));
			// std::cout << "read " << i << " " << npchi << std::endl;
			int parent;
			istrm.read((char*)&parent,sizeof(int));
			// std::cout << "read " << i << " " << parent << std::endl;
			add_rotamer(resn, mychi, npchi, parent);
		}
		std::string endtag;
		istrm >> endtag;
		ALWAYS_ASSERT(endtag=="rot_index_end");
		build_index();
	}

	// assumes primary parent is added firrt among its class
	int
	add_rotamer(
		std::string resname,
		std::vector<float> const & chi,
		int n_proton_chi,
		int parent_key = -1,
		bool is_d_aa = false
	){
		ALWAYS_ASSERT_MSG( n_primary_rotamers_!=0 || parent_key == -1, "must add primary rotamers before children" )
		ALWAYS_ASSERT_MSG( !seen_child_rotamer_ || parent_key != -1, "can't intsert primary rotamer after inserting child" )
		ALWAYS_ASSERT( parent_key == -1 | parent_key < n_primary_rotamers_ );
		Rotamer r;

		rotgen_.get_atoms( resname, chi, r.atoms_, r.hbonders_, r.nheavyatoms, r.donors_, r.acceptors_, d_l_map_ );
		r.resname_ = resname;
		r.chi_ = chi;
		r.n_proton_chi_ = n_proton_chi;
		rotamers_.push_back(r);
		int this_key = rotamers_.size()-1;
		if( parent_key == -1 ){
			n_primary_rotamers_ = this_key+1;
			parent_rotamer_.push_back( this_key );
		} else {
			seen_child_rotamer_ = true;
			parent_rotamer_.push_back( parent_key );
		}
		return this_key;
	}

	void
	build_index()
	{
		std::cout << "n_p_rots: " << n_primary_rotamers_ << std::endl;
		std::cout << "n_rots: " << rotamers_.size() << std::endl;
		std::cout << "n_par_rots: " << parent_rotamer_.size() << std::endl;

		// bounds
		bounds_map_[ rotamers_.front().resname_ ].first = 0;
		for(int i = 1; i < n_primary_rotamers_; ++i){
			if( rotamers_.at(i).resname_ != rotamers_.at(i-1).resname_ ){
				bounds_map_[ resname(i-1) ].second  = i;
				bounds_map_[ resname(i  ) ].first   = i;
			}
		}
		bounds_map_[ rotamers_.back().resname_ ].second = n_primary_rotamers_;
		//utility_exit_with_message("done?");
		// primary rotamers
		// for(int i = 0; i < rotamers_.size(); ++i){
			// if( parent_rotamer_.at(i)==i ) primary_rotamers_.push_back(i);
			// to_primary_index.at(i) = n_primary_rotamers_-1;
		// }

		// assumes all children are contiguous!!!!!!!!!
		child_map_.resize(n_primary_rotamers_,std::make_pair( std::numeric_limits<int>::max(),std::numeric_limits<int>::min() ));
		std::vector<int> child_count(n_primary_rotamers_,0);
		for(int i = 0; i < rotamers_.size(); ++i){
			if( is_primary(i) ) continue; // primary not in list of children

			int ipri = parent_rotamer_.at(i);
			child_count.at(ipri)++;
			child_map_.at(ipri).first  = std::min( child_map_.at(ipri).first , i );
			child_map_.at(ipri).second = std::max( child_map_.at(ipri).second, i );
		}
		//std::cout << "child_map size: " << child_map_.size() << std::endl;
		for(int ipri = 0; ipri < child_map_.size(); ++ipri){
			if( child_count.at(ipri) ){
				++child_map_.at(ipri).second;
			} else {
				child_map_.at(ipri).first = 0;
				child_map_.at(ipri).second = 0;
			}
		}

		for(int i = 0; i < size(); ++i){
			if( chem_index_.have_res( resname(i) ) ) continue;
			for( int ia = 0; ia < natoms(i); ++ia ){
				chem_index_.add_atomdata( resname(i), ia, atom(i,ia).data() );
			}
		}

		protonchi_parent_of_.resize(size(),-1);
		for(int k = 0; k < size(); ++k) protonchi_parent_of_[k] = k;
		for( int irot = 0; irot < size(); ++irot ){
			if( nprotonchi(irot) == 0 ) continue;
			for( int jrot = 0; jrot < irot; ++jrot){
				if( oneletter(irot) != oneletter(jrot) ) continue;
				bool same_noproton_chi = true;
				for( int ichi = 0; ichi < nchi_noproton(irot); ++ichi ){
					same_noproton_chi &= impl::angle_is_close( chi(irot,ichi), chi(jrot,ichi), 0.001f );
				}
				if(same_noproton_chi){
					protonchi_parent_of_.at(irot) = jrot;
					// std::cerr << "pcp " << resname(irot) << " " << irot << " " << jrot;
					// for( int ichi = 0; ichi < nchi(irot); ++ichi){
					// 	std::cerr << " " << chi(irot,ichi) << "/" << chi(jrot,ichi);
					// }
					// std::cerr << std::endl;
					break;
				}
			}
		}

		{
			for( int irot = 0; irot < size(); ++irot ){
				ALWAYS_ASSERT( protonchi_parent_of_.at(irot) >= 0 && protonchi_parent_of_.at(irot) <= irot );
			}
			int uniq_sum = 0;
			for( int i = 0; i < this->size(); ++i ){
				uniq_sum += ( i == protonchi_parent_of_.at(i) );
			}
			std::cerr << "total num protonchi_parent_of_ " << uniq_sum << std::endl;
		}

		structural_parent_of_.resize(size(),-1);
		for(int k = 0; k < size(); ++k) structural_parent_of_[k] = k;
		for( int irot = 0; irot < size(); ++irot ){
			for( int jrot = 0; jrot < irot; ++jrot){
				if( oneletter(irot) != oneletter(jrot) ) continue;
				bool same_chi34 = true;
				for( int ichi = same_struct_start_chi(irot); ichi < nchi_noproton(irot); ++ichi ){
					same_chi34 &= impl::angle_is_close( chi(irot,ichi), chi(jrot,ichi), 0.001f );
				}
				if(same_chi34){
					structural_parent_of_.at(irot) = jrot;
					// std::cerr << "chi34 " << resname(irot) << " " << irot << " " << jrot;
					// for( int ichi = 0; ichi < nchi(irot); ++ichi){
					// 	std::cerr << " " << chi(irot,ichi) << "/" << chi(jrot,ichi);
					// }
					// std::cerr << std::endl;
					break;
				}
			}
		}
		for( int irot = 0; irot < size(); ++irot ){
			if(structural_parent_of_.at(irot) == irot ) structural_parents_.push_back(irot);
		}
		{
			for( int irot = 0; irot < size(); ++irot ){
				ALWAYS_ASSERT( structural_parent_of_.at(irot) >= 0 && structural_parent_of_.at(irot) <= irot );
			}
			int uniq_sum = 0;
			for( int irot = 0; irot < size(); ++irot ){
				uniq_sum += ( irot == structural_parent_of_.at(irot) );
				// if( structural_parent_of_.at(irot)==irot ){
				// 	std::cerr << "chi12_parent: " << irot << " " << resname(irot);
				// 	for( int ichi = 0; ichi < nchi(irot); ++ichi){
				// 		std::cerr << " " << chi(irot,ichi);
				// 	}
				// 	// std::cerr << std::endl;
				// 	int nchild = 0;
				// 	for( int jrot = irot+1; jrot < size(); ++jrot ){
				// 		if( structural_parent_of_.at(jrot) == irot ){
				// 			// std::cerr << "       child: " << jrot << " " << resname(jrot);
				// 			// for( int ichi = 0; ichi < nchi(jrot); ++ichi){
				// 				// std::cerr << " " << chi(jrot,ichi);
				// 			// }
				// 			// std::cerr << std::endl;
				// 			nchild++;
				// 		}
				// 	}
				// 	std::cerr << " nchild: " << nchild << std::endl;
				// }
			}
			std::cerr << "total num structural_parents_ " << structural_parents_.size() << std::endl;
		}

		to_structural_parent_frame_.resize(size(),Xform::Identity());
		for( int irot = 0; irot < size(); ++irot ){
			int isp = structural_parent_of_.at(irot);
			if( isp != irot ){ // if same, leave as identity
				Xform xparent = make_sidechain_stub(isp);
				Xform xthis   = make_sidechain_stub(irot);
				to_structural_parent_frame_[irot] = xparent * xthis.inverse();
			}
		}


		for( int i = 0; i < this->size(); ++i ){
			if( resname(i) == "ALA" ){
				ala_rot_ = i;
				break;
			}
		}

		is_d_.resize( this->size() );
		for ( size_t irot = 0; irot < this->size(); irot++ ) {
			bool is_d =  d_l_map_.count(resname(irot)) != 0;
			is_d_[irot] = is_d;
		}

		sanity_check();
	}

	void
	sanity_check() const
	{
		using std::cerr;
		using std::endl;

		// cerr << "SANITY_CHECK" << endl;
		for( auto const & tmp : bounds_map_ ){
			// cerr << tmp.first << " " << tmp.second.first << " " << tmp.second.second << endl;
			ALWAYS_ASSERT_MSG( tmp.second.first <= tmp.second.second, "RotamerIndex::sanity_check FAIL" );
		}

		for( int irot = 0; irot < size(); ++irot ){
			for( int jrot = 0; jrot < irot; ++jrot ){
				if( resname(irot) != resname(jrot) ) continue;
				bool duplicate_rotamer = true;
				for(int ichi = 0; ichi < nchi(irot); ++ichi){
					duplicate_rotamer &= impl::angle_is_close( chi(irot,ichi), chi(jrot,ichi), 5.0f );
				}
				if( duplicate_rotamer ){
					cerr << "duplicate_rotamer" << endl;
					cerr << irot << " " << resname(irot);
					for(int ichi = 0; ichi < nchi(irot); ++ichi) cerr << " " << chi(irot,ichi); cerr << endl;
					cerr << jrot << " " << resname(jrot);
					for(int ichi = 0; ichi < nchi(jrot); ++ichi) cerr << " " << chi(jrot,ichi); cerr << endl;
				}
				ALWAYS_ASSERT( !duplicate_rotamer )
			}
		}

		for( int irot = 0; irot < size(); ++irot ){
			int ipri = parent_rotamer_.at(irot);
			ALWAYS_ASSERT_MSG( resname(irot) == resname(ipri), "RotamerIndex::sanity_check FAIL" );
			ALWAYS_ASSERT_MSG( rotamers_.at(irot).chi_.size() == rotamers_.at(ipri).chi_.size(), "RotamerIndex::sanity_check FAIL" );
			ALWAYS_ASSERT( rotamers_.at(ipri).n_proton_chi_ <= 1 );
			int nchi = rotamers_.at(ipri).chi_.size() - rotamers_.at(ipri).n_proton_chi_;
			for(int ichi = 0; ichi < nchi; ++ichi){
				float child_chi = chi(irot,ichi);
				float parent_chi = chi(ipri,ichi);
				float chidiff = std::fabs(parent_chi-child_chi);
				chidiff = std::min( chidiff, 360.0f-chidiff );
				// cerr << resname(irot) << " " << ichi << " " << parent_chi << " " << child_chi << " " << chidiff << endl;
				if( ichi==nchi-1 ){
					ALWAYS_ASSERT_MSG( chidiff < 30.0, "parent chi more than 30° from child! "
						+resname(irot)+", "+str(parent_chi)+" vs "+str(child_chi) );
				} else {
					ALWAYS_ASSERT_MSG( chidiff < 20.0, "parent chi more than 20° from child! "
						+resname(irot)+", "+str(parent_chi)+" vs "+str(child_chi) );
				}
			}
		}

		for( int irot = 0; irot < size(); ++irot ){
			if( is_primary(irot) ) continue;
			int ipri = parent_rotamer_.at(irot);
			ALWAYS_ASSERT( child_map_.at(ipri).first <= irot );
			ALWAYS_ASSERT( child_map_.at(ipri).second > irot )
		}

		for( int ipri = 0; ipri < n_primary_rotamers_; ++ipri ){
			for( int ichild = child_map_.at(ipri).first; ichild < child_map_.at(ipri).second; ++ichild ){
				ALWAYS_ASSERT( parent_rotamer_.at(ichild) == ipri )
				ALWAYS_ASSERT( resname(ipri) == resname(ichild) )
				ALWAYS_ASSERT( nchi(ipri) == nchi(ichild) )
				ALWAYS_ASSERT( nprotonchi(ipri) == nprotonchi(ichild) )
				ALWAYS_ASSERT( natoms(ipri) == natoms(ichild) )
				ALWAYS_ASSERT( nheavyatoms(ipri) == nheavyatoms(ichild) )
				for(int ichi = 0; ichi < nchi(ipri); ++ichi){
					ALWAYS_ASSERT( impl::angle_is_close( chi(ipri,ichi), chi(ichild,ichi), 30.0f ) )
				}
			}
		}

		for( int irot = 0; irot < size(); ++irot ){
			for( int ia = 0; ia < rotamers_.at(irot).atoms_.size(); ++ia ){
				AtomData const & ad1( rotamers_.at(irot).atoms_[ia].data() );
				// cerr <<  resname(irot) << endl;
				AtomData const & ad2( chem_index_.atom_data( resname(irot), ia ) );
				ALWAYS_ASSERT( ad1 == ad2 );
			}
		}

		for( int irot = 0; irot < size(); ++irot ){
			int ipcp = protonchi_parent_of_.at(irot);
			ALWAYS_ASSERT( 0 <= protonchi_parent_of_.at(irot) && ipcp < this->size() );
			ALWAYS_ASSERT( resname(ipcp) == resname(irot) )
			ALWAYS_ASSERT( nchi(irot) == nchi(ipcp) )
			ALWAYS_ASSERT( nprotonchi(irot) == nprotonchi(ipcp) )
			ALWAYS_ASSERT( natoms(irot) == natoms(ipcp) )
			ALWAYS_ASSERT( nheavyatoms(irot) == nheavyatoms(ipcp) )
			for(int ichi = 0; ichi < nchi_noproton(irot); ++ichi){
				ALWAYS_ASSERT( impl::angle_is_close( chi(irot,ichi), chi(ipcp,ichi), 0.001f ) )
			}
		}

		for( int irot = 0; irot < size(); ++irot ){
			int isp = structural_parent_of_.at(irot);
			ALWAYS_ASSERT( 0 <= structural_parent_of_.at(irot) && isp < this->size() );
			ALWAYS_ASSERT( resname(isp) == resname(irot) || resname(irot)=="HIS_D" )
			ALWAYS_ASSERT( nchi(irot) == nchi(isp) )
			ALWAYS_ASSERT( nprotonchi(irot) == nprotonchi(isp) )
			ALWAYS_ASSERT( natoms(irot) == natoms(isp) )
			ALWAYS_ASSERT( nheavyatoms(irot) == nheavyatoms(isp) )
			for(int ichi = same_struct_start_chi(irot); ichi < nchi_noproton(irot); ++ichi){
				ALWAYS_ASSERT( impl::angle_is_close( chi(irot,ichi), chi(isp,ichi), 0.001f ) )
			}
		}

		ALWAYS_ASSERT( resname(ala_rot_) == "ALA" );

		for( int irot = 0; irot < size(); ++irot ){
			if( irot < n_primary_rotamers_ ){ ALWAYS_ASSERT(  is_primary(irot) ) }
			else                            { ALWAYS_ASSERT( !is_primary(irot) ) }
		}
	}

	float
	get_max_nbr_radius() const {
		runtime_assert(per_thread_rotamers_.size() > 0);

		float max_radius = 0;
		for ( int irot = 0; irot < size(); irot++ ) {
			max_radius = std::max<float>(max_radius, get_per_thread_rotamer(0, irot)->nbr_radius());
		}
		return max_radius;
	}

	void
	build_per_thread_rotamers( int threads ) {
		per_thread_rotamers_.clear();
		per_thread_lkbrinfo_.clear();
		per_thread_rotamers_.resize(threads);
		per_thread_lkbrinfo_.resize(threads);

		for ( int i = 0; i < per_thread_rotamers_.size(); i++ ) {
			per_thread_rotamers_[i].resize( size() );
			per_thread_lkbrinfo_[i].resize( size() );
		}


		core::chemical::ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");

		for ( int irot = 0; irot < size(); irot++ ) {
			// take care of d vs. l 
			core::chemical::ResidueTypeOP rtypeOP;
			if (d_l_map_.find(resname(irot)) != d_l_map_.end()) {
				core::chemical::ResidueTypeCOP rt_tmp = rts.lock() -> name_map(d_l_map_.find(resname(irot))->second).get_self_ptr();
 				core::chemical::ResidueTypeCOP tmp_dcop = rts.lock() -> get_d_equivalent(rt_tmp);
 				rtypeOP = tmp_dcop -> clone();
 			} else {
 				core::chemical::ResidueTypeCOP rt_tmp = rts.lock() -> name_map(resname(irot)).get_self_ptr();
 				rtypeOP = rt_tmp -> clone();
 			}
			//core::chemical::ResidueType const & rtype = rts.lock()->name_map(resname(irot));
 			core::chemical::ResidueType const & rtype = *rtypeOP;
			core::conformation::ResidueOP rsd = get_residue_at_identity(rtype, chis(irot));
			core::scoring::lkball::LKB_ResidueInfoOP lkbrinfo( new core::scoring::lkball::LKB_ResidueInfo( *rsd ));

			for ( int i = 0; i < per_thread_rotamers_.size(); i++ ) {
				per_thread_rotamers_[i].at(irot) = rsd->clone();
				per_thread_lkbrinfo_[i].at(irot) = std::dynamic_pointer_cast<core::scoring::lkball::LKB_ResidueInfo>(lkbrinfo->clone());
			}
		}
	}

	core::conformation::ResidueOP
	get_per_thread_rotamer( int thread, int irot ) const {
		return per_thread_rotamers_.at(thread).at(irot);
	}
	core::conformation::ResidueOP
	get_per_thread_rotamer_at_identity( int thread, int irot ) const {
		core::conformation::ResidueOP resop = get_per_thread_rotamer( thread, irot );
		residue_to_identity( resop );
		return resop;
	}
	core::scoring::lkball::LKB_ResidueInfoOP
	get_per_thread_lkbrinfo( int thread, int irot ) const {
		return per_thread_lkbrinfo_.at(thread).at(irot);
	}

	std::vector<float> const & chis(int rotnum) const { return rotamers_.at(rotnum).chi_; }
	std::vector<Atom> const & atoms(int rotnum) const { return rotamers_.at(rotnum).atoms_; }

	int ala_rot() const { return ala_rot_; }

	// This function is broken if you add extra residues to the rotamer spec (like the hotspot code does)
	// I don't think you can ever rely on this because of that
	std::pair<int,int>
	index_bounds( std::string const & resname ) const {
		BoundsMap::const_iterator i = bounds_map_.find(resname);
		if( i == bounds_map_.end() ) return std::make_pair(0,0);
		return i->second;
	}

	std::pair<int,int>
	child_bounds_of_primary( int ipri ) const {
		return child_map_.at(ipri);
	}

	Xform
	make_sidechain_stub(int irot) const {
		// std::cerr << "make_sidechain_stub " << irot << " "<< nheavyatoms(irot) << std::endl;
		int const n = nheavyatoms(irot);
		auto p0 = rotamers_.at(irot).atoms_.at(n-1).position();
		auto p1 = rotamers_.at(irot).atoms_.at(n-2).position();
		auto p2 = rotamers_.at(irot).atoms_.at(n-3).position();
		return ::scheme::chemical::make_stub<Xform>(p0, p1, p2);
	}

	//bool is_primary( int irot ) const { return parent_rotamer_.at(irot)==irot; }
	bool is_primary( int irot ) const { return irot < n_primary_rotamers_; }
	bool is_structural_primary( int irot ) const { return structural_parent_of_.at(irot) == irot; }

	void dump_pdb( std::ostream & out, int irot, Xform x=Xform::Identity() ) const {
		// std::cerr << resname(irot) << " " << irot;
		// for(int i = 0; i < rotamers_.at(irot).chi_.size(); ++i) std::cerr << " " << rotamers_.at(irot).chi_.at(i);
		// std::cerr << std::endl;
		out << "MODEL " << resname(irot) << " " << irot << std::endl;
		for( auto const & a : rotamers_.at(irot).atoms_ ){
			auto a2 = a;
			a2.set_position( x * a2.position() );
			std::string s = scheme::io::dump_pdb_atom(a2);
			if( s.size() > 0 && s.size() < 999 ) out << s << std::endl;
		}
		// out << "ENDMDL" << irot << std::endl;
		// out << "MODEL " << resname(irot) << " " << ++rescount << " " << irot << " HBONDERS" << std::endl;
		// for( auto const & h : rotamers_.at(irot).hbonders_ ){
		// 	out << scheme::io::dump_pdb_atom(h.first) << std::endl;
		// 	out << scheme::io::dump_pdb_atom(h.second) << std::endl;
		// }
		for( auto const & hr : rotamers_.at(irot).donors_ ){
			scheme::io::dump_pdb_atom_resname_atomname( out, "DON", "CDON", x*(hr.horb_cen) );
			scheme::io::dump_pdb_atom_resname_atomname( out, "DON", "DDON", x*(hr.horb_cen+hr.direction) );
		}
		for( auto const & hr : rotamers_.at(irot).acceptors_ ){
			scheme::io::dump_pdb_atom_resname_atomname( out, "ACC", "CACC", x*(hr.horb_cen) );
			scheme::io::dump_pdb_atom_resname_atomname( out, "ACC", "DACC", x*(hr.horb_cen+hr.direction) );
		}
		out << "ENDMDL" << std::endl;
	}

	void dump_pdb( std::ostream & out, std::string resn="" ) const {
		std::pair<int,int> b(0,rotamers_.size());
		if( resn.size() ) b = index_bounds(resn);
		int rescount = 0;
		for(int irot = b.first; irot < b.second; ++irot){
			dump_pdb( out, irot );
		}
	}

	void dump_pdb_with_children( std::ostream & out, int iparent ) const {
		ALWAYS_ASSERT( iparent < n_primary_rotamers_ )
		dump_pdb(out, iparent );
		for( int ichild = child_map_.at(iparent).first; ichild < child_map_.at(iparent).second; ++ichild ){
			dump_pdb(out, ichild);
		}
	}

	void dump_pdb_by_structure( std::ostream & out, int isp ) const {
		for( int isc = 0; isc < size(); ++isc ){
			if( structural_parent_of_.at(isc) == isp ){
				dump_pdb( out, isc, to_structural_parent_frame_.at(isc) );
			}
		}
	}

	// template< class Xform >
	void dump_pdb( std::ostream & out, int irot, Xform x, int ires ) const {
		for( int ia = 0; ia < rotamers_.at(irot).atoms_.size(); ++ia ){
			io::dump_pdb_atom( out, x*rotamers_.at(irot).atoms_[ia].position(), rotamers_.at(irot).atoms_[ia].data(), ires );
		}
		for( auto const & h : rotamers_.at(irot).hbonders_ ){
			// io::dump_pdb_atom( out, x*h.first .position(), h.first .data(), ires );
			if( h.second.type() < 0 ){ // is acceptor orbital
				io::dump_pdb_atom( out, x*h.second.position(), h.second.data(), ires );
			}
		}
	}

	uint64_t validation_hash() const {
		uint64_t h = 0;
		for( Rotamer const & rot : rotamers_  ){
			h ^= rot.validation_hash();
		}
		return h;
	}

	int protonchi_parent( int i ) const { return protonchi_parent_of_.at(i); }
	int structural_parent( int i ) const { return structural_parent_of_.at(i); }


	void fill_oneletter_map( std::map<std::string,std::string> & oneletter_map ){
		oneletter_map["ALA"] = 'A';
		oneletter_map["CYS"] = 'C';
		oneletter_map["ASP"] = 'D';
		oneletter_map["GLU"] = 'E';
		oneletter_map["PHE"] = 'F';
		oneletter_map["GLY"] = 'G';
		oneletter_map["HIS"] = 'H';
		oneletter_map["HIS_D"] = 'H';
		oneletter_map["ILE"] = 'I';
		oneletter_map["LYS"] = 'K';
		oneletter_map["LEU"] = 'L';
		oneletter_map["MET"] = 'M';
		oneletter_map["ASN"] = 'N';
		oneletter_map["PRO"] = 'P';
		oneletter_map["GLN"] = 'Q';
		oneletter_map["ARG"] = 'R';
		oneletter_map["SER"] = 'S';
		oneletter_map["THR"] = 'T';
		oneletter_map["VAL"] = 'V';
		oneletter_map["TRP"] = 'W';
		oneletter_map["TYR"] = 'Y';
		oneletter_map["ADE"] = 'a';
		oneletter_map["CYT"] = 'c';
		oneletter_map["GUA"] = 'g';
		oneletter_map["THY"] = 't';
		oneletter_map["RAD"] = 'a';
		oneletter_map["RCY"] = 'c';
		oneletter_map["RGU"] = 'g';
		oneletter_map["URA"] = 'u';
		oneletter_map["H2O"] = 'w';
		oneletter_map["UNP"] = 'z';
		oneletter_map["UNK"] = 'Z';
		oneletter_map["VRT"] = 'X';
		oneletter_map["ala"] = 'A';
		oneletter_map["cys"] = 'C';
		oneletter_map["asp"] = 'D';
		oneletter_map["glu"] = 'E';
		oneletter_map["phe"] = 'F';
		oneletter_map["gly"] = 'G';
		oneletter_map["his"] = 'H';
		oneletter_map["his_d"] = 'H';
		oneletter_map["ile"] = 'I';
		oneletter_map["lys"] = 'K';
		oneletter_map["leu"] = 'L';
		oneletter_map["met"] = 'M';
		oneletter_map["asn"] = 'N';
		oneletter_map["pro"] = 'P';
		oneletter_map["gln"] = 'Q';
		oneletter_map["arg"] = 'R';
		oneletter_map["ser"] = 'S';
		oneletter_map["thr"] = 'T';
		oneletter_map["val"] = 'V';
		oneletter_map["trp"] = 'W';
		oneletter_map["tyr"] = 'Y';
		oneletter_map["ade"] = 'a';
		oneletter_map["cyt"] = 'c';
		oneletter_map["gua"] = 'g';
		oneletter_map["thy"] = 't';
		oneletter_map["rad"] = 'a';
		oneletter_map["rcy"] = 'c';
		oneletter_map["rgu"] = 'g';
		oneletter_map["ura"] = 'u';
		oneletter_map["h2o"] = 'w';
		oneletter_map["unp"] = 'z';
		oneletter_map["unk"] = 'Z';
		oneletter_map["vrt"] = 'X';
		// oneletter_map["DAL"] = 'a';
		// oneletter_map["DCS"] = 'c';
		// oneletter_map["DAS"] = 'd';
		// oneletter_map["DGU"] = 'e';
		// oneletter_map["DPH"] = 'f';
		// oneletter_map["DHI"] = 'h';
		// oneletter_map["DIL"] = 'i';
		// oneletter_map["DLY"] = 'k';
		// oneletter_map["DLE"] = 'l';
		// oneletter_map["DME"] = 'm';
		// oneletter_map["DAN"] = 'n';
		// oneletter_map["DPR"] = 'p';
		// oneletter_map["DGN"] = 'q';
		// oneletter_map["DAR"] = 'r';
		// oneletter_map["DSE"] = 's';
		// oneletter_map["DTH"] = 't';
		// oneletter_map["DVA"] = 'v';
		// oneletter_map["DTR"] = 'w';
		// oneletter_map["DTY"] = 'y';

		oneletter_map["DAL"] = "DAL";
		oneletter_map["DCS"] = "DCS";
		oneletter_map["DAS"] = "DAS";
		oneletter_map["DGU"] = "DGU";
		oneletter_map["DPH"] = "DPH";
		oneletter_map["DHI"] = "DHI";
		oneletter_map["DIL"] = "DIL";
		oneletter_map["DLY"] = "DLY";
		oneletter_map["DLE"] = "DLE";
		oneletter_map["DME"] = "DME";
		oneletter_map["DAN"] = "DAN";
		oneletter_map["DPR"] = "DPR";
		oneletter_map["DGN"] = "DGN";
		oneletter_map["DAR"] = "DAR";
		oneletter_map["DSE"] = "DSE";
		oneletter_map["DTH"] = "DTH";
		oneletter_map["DVA"] = "DVA";
		oneletter_map["DTR"] = "DTR";
		oneletter_map["DTY"] = "DTY";
	}
	void fill_d_l_map( std::map<std::string,std::string> & d_l_map){
		d_l_map.insert(std::make_pair("DAL","ALA"));
		d_l_map.insert(std::make_pair("DCS","CYS"));
		d_l_map.insert(std::make_pair("DAS","ASP"));
		d_l_map.insert(std::make_pair("DGU","GLU"));
		d_l_map.insert(std::make_pair("DPH","PHE"));
		d_l_map.insert(std::make_pair("DHI","HIS"));
		d_l_map.insert(std::make_pair("DIL","ILE"));
		d_l_map.insert(std::make_pair("DLY","LYS"));
		d_l_map.insert(std::make_pair("DLE","LEU"));
		d_l_map.insert(std::make_pair("DME","MET"));
		d_l_map.insert(std::make_pair("DAN","ASN"));
		d_l_map.insert(std::make_pair("DPR","PRO"));
		d_l_map.insert(std::make_pair("DGN","GLN"));
		d_l_map.insert(std::make_pair("DAR","ARG"));
		d_l_map.insert(std::make_pair("DSE","SER"));
		d_l_map.insert(std::make_pair("DTH","THR"));
		d_l_map.insert(std::make_pair("DVA","VAL"));
		d_l_map.insert(std::make_pair("DTR","TRP"));
		d_l_map.insert(std::make_pair("DTY","TYR"));
		//d_l_map.insert(std::make_pair("DHI","HIS_D"));
		//d_l_map.insert(std::make_pair("DGL","GLY"));
	}
	void d_l_map_reverse(std::string l_name, std::string &d_name) {
        for (auto it : this -> d_l_map_) {
        	if (it.second == l_name) {
        		d_name = it.first;
        		break;
        	} else {
        		d_name = l_name;
        	}
        }
	}

	bool operator==(THIS const & other) const {
		return (
				n_primary_rotamers_         == other.n_primary_rotamers_   &&
				seen_child_rotamer_         == other.seen_child_rotamer_   &&
				ala_rot_                    == other.ala_rot_              &&
				chem_index_                 == other.chem_index_           &&
				rotamers_                   == other.rotamers_             &&
				parent_rotamer_             == other.parent_rotamer_       &&
				bounds_map_                 == other.bounds_map_           &&
				child_map_                  == other.child_map_            &&
				protonchi_parent_of_        == other.protonchi_parent_of_  &&
				structural_parents_         == other.structural_parents_   &&
				structural_parent_of_       == other.structural_parent_of_ // &&
				// to_structural_parent_frame_ == other.to_structural_parent_frame_
			);
	}
};


template<class A, class RG, class X>
std::ostream & operator << ( std::ostream & out, RotamerIndex<A,RG,X> const & ridx ){
	out << "RotamerIndex:" << std::endl;
	std::pair<int,int> ib;
	ib=ridx.index_bounds("ALA"); out<<"    ALA "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("CYS"); out<<"    CYS "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("ASP"); out<<"    ASP "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("GLU"); out<<"    GLU "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("PHE"); out<<"    PHE "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("GLY"); out<<"    GLY "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("HIS"); out<<"    HIS "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("ILE"); out<<"    ILE "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("LYS"); out<<"    LYS "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("LEU"); out<<"    LEU "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("MET"); out<<"    MET "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("ASN"); out<<"    ASN "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("PRO"); out<<"    PRO "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("GLN"); out<<"    GLN "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("ARG"); out<<"    ARG "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("SER"); out<<"    SER "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("THR"); out<<"    THR "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("VAL"); out<<"    VAL "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("TRP"); out<<"    TRP "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("TYR"); out<<"    TYR "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;

	ib=ridx.index_bounds("DAL"); out<<"    DAL "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("DCS"); out<<"    DCS "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("DAS"); out<<"    DAS "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("DGU"); out<<"    DGU "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("DPH"); out<<"    DPH "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("DHI"); out<<"    DHI "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("DIL"); out<<"    DIL "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("DLY"); out<<"    DLY "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("DLE"); out<<"    DLE "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("DME"); out<<"    DME "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("DAN"); out<<"    DAN "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("DPR"); out<<"    DPR "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("DGN"); out<<"    DGN "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("DAR"); out<<"    DAR "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("DSE"); out<<"    DSE "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("DTH"); out<<"    DTH "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("DVA"); out<<"    DVA "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("DTR"); out<<"    DTR "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	ib=ridx.index_bounds("DTY"); out<<"    DTY "<<(ib.second-ib.first)<<" "<<ib.first<<"-"<<ib.second-1<<" "<<ridx.nchi(ib.first)<<" "<<ridx.nprotonchi(ib.first)<<std::endl;
	
 	return out;
}





}}

#endif

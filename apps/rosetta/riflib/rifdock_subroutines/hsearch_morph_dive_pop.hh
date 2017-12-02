// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:


#ifndef INCLUDED_riflib_rifdock_subroutines_hsearch_morph_dive_pop_hh
#define INCLUDED_riflib_rifdock_subroutines_hsearch_morph_dive_pop_hh


#include <riflib/types.hh>
#include <riflib/rifdock_typedefs.hh>
#include <riflib/rifdock_subroutines/util.hh>
#include <riflib/rifdock_subroutines/meta.hh>
#include <riflib/rifdock_subroutines/hsearch_original.hh>


using ::scheme::make_shared;
using ::scheme::shared_ptr;

typedef int32_t intRot;


template<class DirectorBase, class ScaffoldProvider >
bool
do_an_hsearch(uint64_t start_resl, 
    std::vector< std::vector< tmplSearchPoint<_DirectorBigIndex<DirectorBase>> > > & samples, 
    HsearchData<DirectorBase, ScaffoldProvider > & d) {


    using namespace core::scoring;
        using std::cout;
        using std::endl;
        using namespace devel::scheme;
        typedef numeric::xyzVector<core::Real> Vec;
        typedef numeric::xyzMatrix<core::Real> Mat;
        // typedef numeric::xyzTransform<core::Real> Xform;
        using ObjexxFCL::format::F;
        using ObjexxFCL::format::I;
        using devel::scheme::print_header;
        using ::devel::scheme::RotamerIndex;

    typedef ::scheme::util::SimpleArray<3,float> F3;
    typedef ::scheme::util::SimpleArray<3,int> I3;

    typedef _SearchPointWithRots<DirectorBase> SearchPointWithRots;

    typedef _DirectorBigIndex<DirectorBase> DirectorIndex;
    typedef tmplSearchPoint<DirectorIndex> SearchPoint;

    typedef typename ScaffoldProvider::ScaffoldIndex ScaffoldIndex;



using ::scheme::scaffold::BOGUS_INDEX;
using ::scheme::scaffold::TreeIndex;
using ::scheme::scaffold::TreeLimits;




    bool search_failed = false;
    {

        for( int this_stage = 0; this_stage < samples.size(); ++this_stage )
        {
            int iresl = this_stage + start_resl;
            cout << "HSearsh stage " << iresl+1 << " resl " << F(5,2,d.RESLS[this_stage]) << " begin threaded sampling, " << KMGT(samples[this_stage].size()) << " samples: ";
            int64_t const out_interval = samples[this_stage].size()/50;
            std::exception_ptr exception = nullptr;
            std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
            start = std::chrono::high_resolution_clock::now();
            d.total_search_effort += samples[this_stage].size();

            #ifdef USE_OPENMP
            #pragma omp parallel for schedule(dynamic,64)
            #endif
            for( int64_t i = 0; i < samples[this_stage].size(); ++i ){
                if( exception ) continue;
                try {
                    if( i%out_interval==0 ){ cout << '*'; cout.flush(); }
                    DirectorIndex const isamp = samples[this_stage][i].index;

                    ScenePtr tscene( d.scene_pt[omp_get_thread_num()] );
                    d.director->set_scene( isamp, iresl, *tscene );

                    // the real rif score!!!!!!
                    samples[this_stage][i].score = d.objectives[iresl]->score( *tscene );// + tot_sym_score;


                } catch( std::exception const & ex ) {
                    #ifdef USE_OPENMP
                    #pragma omp critical
                    #endif
                    exception = std::current_exception();
                }
            }
            if( exception ) std::rethrow_exception(exception);
            end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed_seconds_rif = end-start;
            float rate = (double)samples[this_stage].size()/ elapsed_seconds_rif.count()/omp_max_threads();
            cout << endl;// << "done threaded sampling, partitioning data..." << endl;

            SearchPoint max_pt, min_pt;
            int64_t len = samples[this_stage].size();
            if( samples[this_stage].size() > d.opt.beam_size/d.opt.DIMPOW2 ){
                __gnu_parallel::nth_element( samples[this_stage].begin(), samples[this_stage].begin()+d.opt.beam_size/d.opt.DIMPOW2, samples[this_stage].end() );
                len = d.opt.beam_size/d.opt.DIMPOW2;
                min_pt = *__gnu_parallel::min_element( samples[this_stage].begin(), samples[this_stage].begin()+len );
                max_pt = *(samples[this_stage].begin()+d.opt.beam_size/d.opt.DIMPOW2);
            } else {
                min_pt = *__gnu_parallel::min_element( samples[this_stage].begin(), samples[this_stage].end() );
                max_pt = *__gnu_parallel::max_element( samples[this_stage].begin(), samples[this_stage].end() );
            }

            cout << "HSearsh stage " << iresl+1 << " complete, resl. " << F(7,3,d.RESLS[iresl]) << ", "
                  << " " << KMGT(samples[this_stage].size()) << ", promote: " << F(9,6,min_pt.score) << " to "
                  << F(9,6, std::min(d.opt.global_score_cut,max_pt.score)) << " rate " << KMGT(rate) << "/s/t " << std::endl;

            // cout << "Answer: " << ( answer_exists ? "exists" : "doesn't exist" ) << std::endl;

            if( this_stage+1 == samples.size() ) break;

            for( int64_t i = 0; i < len; ++i ){
                uint64_t isamp0 = ::scheme::kinematics::bigindex_nest_part(samples[this_stage][i].index);
                if( samples[this_stage][i].score >= d.opt.global_score_cut ) continue;
                if( iresl == 0 ) ++d.non0_space_size;
                for( uint64_t j = 0; j < d.opt.DIMPOW2; ++j ){
                    uint64_t isamp = isamp0 * d.opt.DIMPOW2 + j;
                    samples[this_stage+1].push_back( SearchPoint(DirectorIndex(isamp, ::scheme::kinematics::bigindex_scaffold_index(samples[this_stage][i].index))) );
                }
            }
            if( 0 == samples[this_stage+1].size() ){
                search_failed = true;
                std::cout << "search fail, no valid samples!" << std::endl;
                break;
            }
            samples[this_stage].clear();

        }
        if( search_failed ) return false;
        std::cout << "full sort of final samples" << std::endl;
        __gnu_parallel::sort( samples.back().begin(), samples.back().end() );
    }
    if( search_failed ) return false;

    return true;


}


// template<__DirectorBase>
// using HsearchFunctionType = typedef


template<class DirectorBase, class ScaffoldProvider >
bool
hsearch_morph_dive_pop(
    shared_ptr<std::vector< _SearchPointWithRots<DirectorBase> > > & hsearch_results_p,
    HsearchData<DirectorBase, ScaffoldProvider > & d) {


    using namespace core::scoring;
        using std::cout;
        using std::endl;
        using namespace devel::scheme;
        typedef numeric::xyzVector<core::Real> Vec;
        typedef numeric::xyzMatrix<core::Real> Mat;
        // typedef numeric::xyzTransform<core::Real> Xform;
        using ObjexxFCL::format::F;
        using ObjexxFCL::format::I;
        using devel::scheme::print_header;
        using ::devel::scheme::RotamerIndex;

    typedef ::scheme::util::SimpleArray<3,float> F3;
    typedef ::scheme::util::SimpleArray<3,int> I3;

    typedef _SearchPointWithRots<DirectorBase> SearchPointWithRots;

    typedef _DirectorBigIndex<DirectorBase> DirectorIndex;
    typedef tmplSearchPoint<DirectorIndex> SearchPoint;

    typedef typename ScaffoldProvider::ScaffoldIndex ScaffoldIndex;


using ::scheme::scaffold::BOGUS_INDEX;
using ::scheme::scaffold::TreeIndex;
using ::scheme::scaffold::TreeLimits;

    shared_ptr<MorphingScaffoldProvider> morph_provider = std::dynamic_pointer_cast<MorphingScaffoldProvider>(d.scaffold_provider);

    morph_provider->test_make_children( TreeIndex(0, 0) );

    TreeLimits limits = morph_provider->get_scaffold_index_limits();
    uint64_t num_scaffolds = limits[1];

    for ( uint64_t scaffno = 0; scaffno < num_scaffolds; scaffno++ ) {
        TreeIndex ti(1, scaffno);
        ScaffoldDataCacheOP sdc = morph_provider->get_data_cache_slow(ti);
        sdc->setup_onebody_tables( d.rot_index_p, d.opt);
    }



    std::vector< std::vector< SearchPoint > > samples( d.RESLS.size() );
    samples[0].resize( ::scheme::kinematics::bigindex_nest_part(d.director->size(0))*num_scaffolds );

    uint64_t index_count = 0;
    for ( uint64_t scaffno = 0; scaffno < num_scaffolds; scaffno++ ) {
        for( uint64_t i = 0; i < ::scheme::kinematics::bigindex_nest_part(d.director->size(0)); ++i ) {
            samples[0][index_count++] = SearchPoint( DirectorIndex( i, TreeIndex(1, scaffno)) );
        }
    }

    BOOST_FOREACH( ScenePtr & s, d.scene_pt ) s = d.scene_minimal->clone_specific_deep(std::vector<uint64_t> {1});



    bool success = do_an_hsearch( 0, samples, d );

    if ( ! success ) return false;


    std::cout << "total non-0 space size was approx " << float(d.non0_space_size)*1024.0*1024.0*1024.0 << " grid points" << std::endl;
    std::cout << "total search effort " << KMGT(d.total_search_effort) << std::endl;


    hsearch_results_p = make_shared<std::vector< SearchPointWithRots >>();
    std::vector< SearchPointWithRots > & hsearch_results = *hsearch_results_p;


    hsearch_results.resize( samples.back().size() );
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,1024)
    #endif
    for( int ipack = 0; ipack < hsearch_results.size(); ++ipack ){
        hsearch_results[ipack].score = samples.back()[ipack].score;
        hsearch_results[ipack].index = samples.back()[ipack].index;
    }



    return true;



}




#endif
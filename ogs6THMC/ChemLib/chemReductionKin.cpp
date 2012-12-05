/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file chemReductionKin.cpp
 *
 * Created on 2012-08-23 by Haibing Shao
 */

#include "chemReductionKin.h"
#include "logog.hpp"

namespace ogsChem
{

chemReductionKin::chemReductionKin(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp, 
	                               std::vector<ogsChem::chemReactionKin*> & list_kin_reactions)
    : _list_kin_reactions(list_kin_reactions)
{
	// by default, the class is not yet initialized
	isInitialized = false; 

	if ( map_chemComp.size() > 0 && _list_kin_reactions.size() > 0 )
	{   // if there are reactions. 
		// cout how many mobile and how many immobile components
		countComp(map_chemComp); 
		// make stoichiometric matrix
		buildStoi(map_chemComp, list_kin_reactions); 
		// calculate the reduction parameters	
		update_reductionScheme(); 
		// flip the initialization flag
		isInitialized = true; 
	}
}

chemReductionKin::~chemReductionKin()
{}

void chemReductionKin::buildStoi(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp, 
	                             std::vector<ogsChem::chemReactionKin*> & list_kin_reactions)
{
	size_t i,j, tmp_idx; 
	double tmp_stoi; 
	std::string tmp_str;
	BaseLib::OrderedMap<std::string, ogsChem::ChemComp*>::iterator tmp_Comp;
	
	// obtain the size info
	_I = map_chemComp.size(); 
	_J = list_kin_reactions.size(); 
	// creat the memory for Stoi matrix
	_matStoi = LocalMatrix::Zero(_I, _J); 
	
	// based on the reactions, fill in the stoi matrix
	// loop over all reactions
	for ( j=0; j < list_kin_reactions.size(); j++ )
	{	// for each reaction
		// find each participating components 
		for ( i=0; i < list_kin_reactions[j]->get_vecCompNames().size(); i++ ){
			tmp_str  = list_kin_reactions[j]->get_vecCompNames()[i]; 
			tmp_Comp = map_chemComp.find(tmp_str);
			tmp_idx  = tmp_Comp->second->getIndex(); 
			if ( list_kin_reactions[j]->get_vecStoi().size() > 2 )
			{
				// normal reactions
				tmp_stoi = list_kin_reactions[j]->get_vecStoi()[i];
				// and put them into Stoi matrix
				_matStoi(tmp_idx,j) = tmp_stoi; 
			}
			else  // this is a basis component
			{
			    _matStoi(tmp_idx,j) = 1.0; 
			}

		}  // end of for i
	}  // end of for j

#ifdef _DEBUG
	// debugging--------------------------
	std::cout << "Stoichiometric Matrix S: " << std::endl; 
	std::cout << _matStoi << std::endl;
	// end of debugging-------------------
#endif
}

void chemReductionKin::update_reductionScheme(void)
{
	// devide mobile and immobile parts
	// S1 =S(1:I,:);
	_matS_1 = _matStoi.topRows( _I_mob );
	// S2 =S(I+1:length(C),:);
	_matS_2 = _matStoi.bottomRows( _I_min ); 

	// reveal the rank of S1 and S2
	Eigen::FullPivLU<LocalMatrix> lu_decomp_S1(_matS_1); 
	_mat_s_1 = lu_decomp_S1.image(_matS_1); 
	Eigen::FullPivLU<LocalMatrix> lu_decomp_S2(_matS_2); 
	_mat_s_2 = lu_decomp_S2.image(_matS_2); 

#ifdef _DEBUG
    std::cout << "_mat_S_1: "    << std::endl; 
	std::cout << _matS_1 << std::endl;
	std::cout << "_mat_S_2: "    << std::endl; 
	std::cout << _matS_2 << std::endl;
	std::cout << "_mat_s_1: "    << std::endl; 
	std::cout << _mat_s_1 << std::endl;
	std::cout << "_mat_s_2: "    << std::endl; 
	std::cout << _mat_s_2 << std::endl;
#endif

	// Calculate the s1T = S_i^T matrix consisting of a max set of linearly
	// independent columns that are orthogonal to each column of S_i^*.
	// s1T = orthcomp(s1);
	_matS_1_ast = orthcomp( _mat_s_1 );
    // s2T = orthcomp(s2);
	_matS_2_ast = orthcomp( _mat_s_2 );

#ifdef _DEBUG
	std::cout << "s_1*: "    << std::endl; 
	std::cout << _matS_1_ast << std::endl;
	std::cout << "s_2*: "    << std::endl; 
	std::cout << _matS_2_ast << std::endl;
#endif

	// now store the corresponding size
	this->_n_eta_mob   = _matS_1_ast.cols();
	this->_n_eta_immob = _matS_2_ast.cols();
	this->_n_xi_mob    = _mat_s_1.cols(); 
	this->_n_xi_immob  = _mat_s_2.cols();
	this->_n_eta       = _n_eta_mob + _n_eta_immob; 
	this->_n_xi        = _n_xi_mob  + _n_xi_immob;

	// now calculated A1 and A2
	_matA1 = ( _mat_s_1.transpose() * _mat_s_1 ).fullPivHouseholderQr().solve(_mat_s_1.transpose()) * _matS_1 ; 
	_matA2 = ( _mat_s_2.transpose() * _mat_s_2 ).fullPivHouseholderQr().solve(_mat_s_2.transpose()) * _matS_2 ; 

#ifdef _DEBUG
	std::cout << "A_1: "    << std::endl; 
	std::cout << _matA1 << std::endl;
	std::cout << "A_2: "    << std::endl; 
	std::cout << _matA2 << std::endl;
#endif

	_mat_c_mob_2_eta_mob     = ( _matS_1_ast.transpose() * _matS_1_ast ).fullPivHouseholderQr().solve(_matS_1_ast.transpose()); 
    _mat_c_immob_2_eta_immob = ( _matS_2_ast.transpose() * _matS_2_ast ).fullPivHouseholderQr().solve(_matS_2_ast.transpose()); 
	_mat_c_mob_2_xi_mob      = ( _mat_s_1.transpose() * _mat_s_1 ).fullPivHouseholderQr().solve(_mat_s_1.transpose());
	_mat_c_immob_2_xi_immob  = ( _mat_s_2.transpose() * _mat_s_2 ).fullPivHouseholderQr().solve(_mat_s_2.transpose()); 

#ifdef _DEBUG
	std::cout << "_mat_c_mob_2_eta_mob: "    << std::endl; 
	std::cout << _mat_c_mob_2_eta_mob << std::endl;
	std::cout << "_mat_c_immob_2_eta_immob: "    << std::endl; 
	std::cout << _mat_c_immob_2_eta_immob << std::endl;
	std::cout << "_mat_c_mob_2_xi_mob: "    << std::endl; 
	std::cout << _mat_c_mob_2_xi_mob << std::endl;
	std::cout << "_mat_c_immob_2_xi_immob: "    << std::endl; 
	std::cout << _mat_c_immob_2_xi_immob << std::endl;
#endif
}

LocalMatrix chemReductionKin::orthcomp( LocalMatrix & inMat )
{
	// initialize it so that they have the same dimension
    LocalMatrix outMat; 

	Eigen::FullPivLU<LocalMatrix> lu_decomp(inMat.transpose());

	outMat = lu_decomp.kernel(); 

	// notice that if the kernel returns a matrix with dimension zero, 
	// then the returned matrix will be a column-vector filled with zeros
	// therefore we do a saftety check here, and set the number of columns to zero
	if ( outMat.cols() == 1 && outMat.col(0).norm() == 0.0 )
		outMat = LocalMatrix::Zero(inMat.rows(), 0);

	return outMat;
}

void chemReductionKin::countComp(BaseLib::OrderedMap<std::string, ogsChem::ChemComp*> & map_chemComp)
{
	_I_mob = 0; 
	_I_sorp= 0; 
	_I_min = 0; 

	BaseLib::OrderedMap<std::string, ogsChem::ChemComp*>::iterator it; 
	for( it = map_chemComp.begin(); it != map_chemComp.end(); it++ )
	{
		switch ( it->second->getMobility() )
		{
		case ogsChem::MOBILE: 
			_I_mob++; 
			break;
		case ogsChem::SORPTION: 
			_I_sorp++;
			break;
		case ogsChem::MINERAL: 
			_I_min++;
			break;
		default:
			_I_min++;
			break; 
		}
	}
}


void chemReductionKin::Conc2EtaXi(ogsChem::LocalVector &local_conc, 
	                              ogsChem::LocalVector &local_eta_mob, 
								  ogsChem::LocalVector &local_eta_immob, 
								  ogsChem::LocalVector &local_xi_mob, 
								  ogsChem::LocalVector &local_xi_immob )
{
	// declare local temp variable
	ogsChem::LocalVector local_c_mob, local_c_immob; 

	// divide c1 and c2
	local_c_mob   = local_conc.topRows(    this->_I_mob );
	local_c_immob = local_conc.bottomRows( this->_I_sorp + this->_I_min ); 

	// convert eta_mob and xi_mob
	local_eta_mob   = _mat_c_mob_2_eta_mob * local_c_mob; 
	local_xi_mob    = _mat_c_mob_2_xi_mob  * local_c_mob; 

	// convert eta_immob and xi_immob
	local_eta_immob = _mat_c_immob_2_eta_immob * local_c_immob;  
    local_xi_immob  = _mat_c_immob_2_xi_immob  * local_c_immob; 


}

void chemReductionKin::EtaXi2Conc(ogsChem::LocalVector &local_eta_mob, 
	                              ogsChem::LocalVector &local_eta_immob, 
								  ogsChem::LocalVector &local_xi_mob, 
								  ogsChem::LocalVector &local_xi_immob,
								  ogsChem::LocalVector &local_conc )
{
	// declare local temp variable
	ogsChem::LocalVector local_c_mob, local_c_immob; 

	local_c_mob   = _mat_s_1 * local_xi_mob   + _matS_1_ast * local_eta_mob; 
	local_c_immob = _mat_s_2 * local_xi_immob + _matS_2_ast * local_eta_immob; 

	local_conc.topRows( this->_I_mob ) = local_c_mob; 
	local_conc.bottomRows( this->_I_sorp + this->_I_min ) = local_c_immob; 

}

void chemReductionKin::Calc_Xi_mob_Rate(ogsChem::LocalVector &local_eta_mob, 
	                                ogsChem::LocalVector &local_eta_immob, 
									ogsChem::LocalVector &local_xi_mob,
									ogsChem::LocalVector &local_xi_immob, 
									ogsChem::LocalVector &xi_mob_rate     )
{
	size_t i; 

	// the size of vec_rates is equal to the number of kinetic equations
	ogsChem::LocalVector vec_rates = ogsChem::LocalVector::Zero(_J); 
	// the local temp concentration vector
	ogsChem::LocalVector vec_conc = ogsChem::LocalVector::Zero(_I); 

	// first convert these eta and xi to concentrations
	EtaXi2Conc( local_eta_mob, 
		        local_eta_immob, 
				local_xi_mob, 
				local_xi_immob, 
				vec_conc ); 

	// then calculate the rates and fill them in the rate vector
	for ( i=0; i < _J; i++ )
	{
		// get to the particular kin equation and calculate its rate
		this->_list_kin_reactions[i]->calcReactionRate( vec_conc ); 
		vec_rates(i) = this->_list_kin_reactions[i]->getRate(); 
	}

	// multiply the rate vector with the A matrix to get rate for xi_mob and xi_immob
	xi_mob_rate   = _matA1 * vec_rates; 
}


void chemReductionKin::Calc_Xi_immob_Rate(ogsChem::LocalVector &local_eta_mob, 
	                                      ogsChem::LocalVector &local_eta_immob, 
									      ogsChem::LocalVector &local_xi_mob,
									      ogsChem::LocalVector &local_xi_immob, 
								          ogsChem::LocalVector &xi_immob_rate     )
{
	size_t i; 

	// the size of vec_rates is equal to the number of kinetic equations
	ogsChem::LocalVector vec_rates = ogsChem::LocalVector::Zero(_J); 
	// the local temp concentration vector
	ogsChem::LocalVector vec_conc = ogsChem::LocalVector::Zero(_I); 

	// first convert these eta and xi to concentrations
	EtaXi2Conc( local_eta_mob, 
		        local_eta_immob, 
				local_xi_mob, 
				local_xi_immob, 
				vec_conc ); 

	// then calculate the rates and fill them in the rate vector
	for ( i=0; i < _J; i++ )
	{
		// get to the particular kin equation and calculate its rate
		this->_list_kin_reactions[i]->calcReactionRate( vec_conc ); 
		vec_rates(i) = this->_list_kin_reactions[i]->getRate(); 
	}

	// multiply the rate vector with the A matrix to get rate for xi_mob and xi_immob
	xi_immob_rate = _matA2 * vec_rates; 
}


void chemReductionKin::Calc_Xi_Rate(ogsChem::LocalVector &local_eta_mob, 
	                                ogsChem::LocalVector &local_eta_immob, 
									ogsChem::LocalVector &local_xi_mob,
									ogsChem::LocalVector &local_xi_immob, 
									ogsChem::LocalVector &xi_mob_rate, 
									ogsChem::LocalVector &xi_immob_rate)
{
	size_t i; 

	// the size of vec_rates is equal to the number of kinetic equations
	ogsChem::LocalVector vec_rates = ogsChem::LocalVector::Zero(_J); 
	// the local temp concentration vector
	ogsChem::LocalVector vec_conc = ogsChem::LocalVector::Zero(_I); 

	// first convert these eta and xi to concentrations
	EtaXi2Conc( local_eta_mob, 
		        local_eta_immob, 
				local_xi_mob, 
				local_xi_immob, 
				vec_conc ); 

	// then calculate the rates and fill them in the rate vector
	for ( i=0; i < _J; i++ )
	{
		// get to the particular kin equation and calculate its rate
		this->_list_kin_reactions[i]->calcReactionRate( vec_conc ); 
		vec_rates(i) = this->_list_kin_reactions[i]->getRate(); 
	}

	// multiply the rate vector with the A matrix to get rate for xi_mob and xi_immob
	xi_mob_rate   = _matA1 * vec_rates; 
	xi_immob_rate = _matA2 * vec_rates; 
}


}  // end of namespace
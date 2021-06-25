import itertools
import numpy as np
import math
import random
import matplotlib.pyplot as plt
import time
import re
from scipy import optimize
from functools import reduce
from operator import *


# Theta_eqs class

class theta_eqs:
    
    PT_mode='scalar'
   
    
    def __init__(self, J, max_PT_order, spin_amount, input_couplings, input_generators, debugging_mode='False'):
        
        stopwatch=time.time()
        
        self.debugging_mode=debugging_mode
        
        self.J=J
        self.max_PT_order=max_PT_order
        self.spin_amount=spin_amount        
        self.computational_states=[list(computational_state) for computational_state in (itertools.product(*[[0,1] for spin in range(spin_amount)]))]
            
        self.couplings=[]

        self.coupling_strengths=[]

        for coupling in input_couplings:

            Tobias_coupling_string=coupling[0]

            coupling_strength=coupling[1]

            assert(coupling_strength==J)

            processed_coupling_string=[[Tobias_coupling_string[2*symbol_id], eval(Tobias_coupling_string[2*symbol_id+1])]
                  for symbol_id in range(len(Tobias_coupling_string)//2)]

            Yaroslav_coupling_string=['I' for spin in range(spin_amount)]

            for sub in processed_coupling_string:

                Yaroslav_coupling_string[sub[1]-1]=sub[0]

            self.couplings+=[Yaroslav_coupling_string]

            self.coupling_strengths+=[coupling_strength]


        self.unitary_generators=[]

        for Tobias_generator_string in input_generators:


            processed_generator_string=[[Tobias_generator_string[2*symbol_id], eval(Tobias_generator_string[2*symbol_id+1])]
                  for symbol_id in range(len(Tobias_generator_string)//2)]

            Yaroslav_generator_string=['I' for spin in range(spin_amount)]

            for sub in processed_generator_string:

                Yaroslav_generator_string[sub[1]-1]=sub[0]

            self.unitary_generators+=[Yaroslav_generator_string]


        

        self.couplings_amount=len(self.couplings)

        self.number_of_thetas=len(self.unitary_generators)

        self.generator_to_theta_dictionary=[[i] for i in range(len(input_generators))]

        assert tuple(sorted([theta for generator_thetas in self.generator_to_theta_dictionary for theta in generator_thetas])
                    ) == tuple(range(self.number_of_thetas))

        assert len(self.generator_to_theta_dictionary)==self.number_of_thetas
        
        self.log(f'theta_eqs object created in {time.time()-stopwatch} seconds',True)
        

# Debugging section

    def log(self, phrase, display_outside_debugging_mode=False):
        '''
        phrase - fstring describing the message
        display_outside_debugging_mode - boolean, explaining whether this line should be printed outside debugging mode
        '''
        if ( (self.debugging_mode==True) | (display_outside_debugging_mode) ):
            print(phrase)
            print()
    
    
    ### k vectors, couplings and generators

    def list_of_Ks_from_PT_order(self, local_vector_PT_mode=False):

    #     to-be internal function:
    # test by running it directly or as a part of a larger procedure
        
        
#        local_PT_mode may be different from self.PT_mode: most notably, vector PT_mode is needed to produce good segregated PT series
        
        if local_vector_PT_mode:
            local_PT_mode='vector'
        else:
            local_PT_mode=self.PT_mode

        if local_PT_mode=='vector':

            list_of_Ks=[]

            for mod_K in range(1, self.max_PT_order+1):
                self.log(f'list_of_Ks addition={list(partitioning(mod_K))}')
                list_of_Ks+=list(partitioning(mod_K))

            list_of_Ks=[fill_the_list(K, self.couplings_amount) for K in list_of_Ks if fill_the_list(K,self.couplings_amount)!=None]

            '''
            including all permutations of the Ks:
            '''

            list_of_Ks=[list(K) for K in list(set(reduce(lambda x, y: x+y, [list(itertools.permutations(K_lexic)) for K_lexic in list_of_Ks])))]

            list_of_Ks=sorted(list_of_Ks, key = lambda K: sum(K))

        elif local_PT_mode=='scalar':

            list_of_Ks=[[PT_order] for PT_order in range(1,self.max_PT_order+1)]

        else:

            raise Exception('Unknown PT_mode!')

        return(list_of_Ks)
    
    def theta_to_generator(self, theta_label):
#     to-be class method:
# test by running it directly or as a part of a larger procedure
    
    
        for thetas_per_generator in self.generator_to_theta_dictionary:

            if theta_label in thetas_per_generator:

                return(self.generator_to_theta_dictionary.index(thetas_per_generator))

        raise Exception('Wrong input: nonexistent theta!')
    
    def pauli_strings_from_generator_list(self, generator_list):

    #     to-be class method:
    # test by running it directly or as a part of a larger procedure

        return([self.unitary_generators[generator_label] for generator_label in generator_list])



    def pauli_strings_from_couplings_list(self, coupling_labels):
    
    # to-be class method:
    # test by running it directly or as a part of a larger procedure

        '''
        Returns a list of pauli strings corresponding to a list of coupling labels

        uses global variable 'couplings'
        '''

        self.log(f'pauli_strings_from_couplings_list: \n The coupling labels are: \n {coupling_labels}; \n The couplings are: \n {self.couplings}')

        return([self.couplings[coupling_label] for coupling_label in coupling_labels])
    
    def function1(self):
        
        return(self.PT_mode)
    
    def function2(self):
        return(self.function1())
    
    ## PT series

    ### Dyson elementaries

    def s_of_k(self, K):

    #     to-be class method:
    # test by running it directly or as a part of a larger procedure

        self.log(f's_of_k:\n K={K},\n k_to_couplings(K) = {k_to_couplings(K)}')

        return(threaded_pauli_strings_action(self.pauli_strings_from_couplings_list(k_to_couplings(K)), [[0 for spin in range(self.spin_amount)],1])[0])
    
    def PT_cutoff(self, C_series):
    
    # to-be class method:
    # test by running it directly or as a part of a larger procedure

        return([C for C in C_series if np.sum(C[0])<=self.max_PT_order])
    
    
    ### Dyson calculus: unnormalized C-series
    
    def K_trivial_state_check(self, K):
    
# remove PT_mode from the arguments (i.e. turn it into a global variable)
# move it to the class body
# test by running it directly or as a part of a larger procedure

        return( reduce(lambda A, B: A and B, [s==0 for s in self.s_of_k(K)] ) and not reduce(lambda A, B: A and B, [Ki==0 for Ki in K] ))
    
    def unnormalized_C_dictionary(self):

    #     to-be class method:
    # move it to the class body
    # test by running it directly or as a part of a larger procedure

        list_of_Ks=self.list_of_Ks_from_PT_order(local_vector_PT_mode=True)

        the_C_list=[[np.array([0 for coupling in range(self.couplings_amount)]), 1]]
        the_C_dictionary={str(a_C[0].tolist()): a_C[1] for a_C in the_C_list}

        delta_betas=[np.array(delta_beta) for delta_beta in np.eye(self.couplings_amount, dtype='int').tolist()]

        for K_as_a_list in [K for K in list_of_Ks if ((not self.K_trivial_state_check(K)) and (sum(K)<=self.max_PT_order))]:
            
            
            self.log(f'unnormalized_C_dictionary iteration:')

            K=np.array(K_as_a_list)

            self.log(f'K={K}, \n delta_betas[0]={delta_betas[0]}')

            K_betas=[K-delta_betas[beta] for beta in range(self.couplings_amount) 
                     if K[beta]>0 and not self.K_trivial_state_check(K-delta_betas[beta])]

            C_betas=[the_C_dictionary[str(K_beta.tolist())] for K_beta in K_betas]

            self.log(f'K_betas={K_betas}')

            k_primes=[np.array(k_prime) for k_prime in list(itertools.product(*[list(range(Ki+1)) for Ki in K]))]

            k_primes=[k_prime for k_prime in k_primes if (self.K_trivial_state_check(k_prime))]

            k_primes_betas=[[k_prime-delta_betas[beta] for beta in range(self.couplings_amount) if k_prime[beta]>0] for k_prime in k_primes]


            self.log(f'k_primes={k_primes}')
            self.log(f'k_primes_betas={k_primes_betas}')

            C_k_prime_beta_sums=[sum([the_C_dictionary[str(k_prime_beta.tolist())] if str(k_prime_beta.tolist()) 
                                      in the_C_dictionary else 0 for k_prime_beta in k_primes_betas[k_prime_index] ]) 
                                 for k_prime_index in range(len(k_primes))]

            C_k_minus_k_primes=[the_C_dictionary[str((K-k_prime).tolist())] if str((K-k_prime).tolist()) 
                                in the_C_dictionary else 0 for k_prime in k_primes]


            k_prime_terms=[C_k_prime_beta_sums[k_prime_iterator]*C_k_minus_k_primes[k_prime_iterator] 
                           for k_prime_iterator in range(len(k_primes))]



            self.log(f'k_primes_terms={k_prime_terms}')

            the_C_list+=[[K, (sum(C_betas)-sum(k_prime_terms))/(E0_of_s([0 for spin in range(self.spin_amount)])

                                                                -E0_of_s(self.s_of_k( K)))]]

            self.log(f'the_C_list addition = {[K, (sum(C_betas)-sum(k_prime_terms))/(E0_of_s([0,0,0,0])-E0_of_s(self.s_of_k(K)))]}')



            the_C_dictionary={str(a_C[0].tolist()): a_C[1] for a_C in the_C_list}

            self.log(f'current the_C_dictionary = {the_C_dictionary}')


        return(the_C_dictionary)


### Normalizing, equation-adapting C-series
    
    
    
    def C_series_to_Z(self, C_series):


    #     to-be class method:
    # test by running it directly or as a part of a larger procedure

        C_series_by_comp_state={str(s): [C for C in C_series if self.s_of_k(C[0])==s] for s in self.computational_states}

        self.log(f'C_series_by_comp_state\n')
        for comp_state in C_series_by_comp_state:
            self.log(f'{comp_state}={C_series_by_comp_state[comp_state]}')

        pre_Z=[C_series_mult(C_series_by_comp_state[str(s)], C_series_by_comp_state[str(s)]) for s in self.computational_states]

        for pre_Z_element in pre_Z:
            self.log(f'pre_Z_element={pre_Z_element}\n')

        Z=reduce(C_series_add, pre_Z)

        self.log(f'Z={Z}')

        return(Z)


    def N_from_Z(self, Z):
        
    #     to-be class method:
    # test by running it directly or as a part of a larger procedure

        X=[C_Z for C_Z in Z if np.sum(C_Z[0])!=0]

        alpha=-1/2

        N=[[np.array([0 for coupling in range(self.couplings_amount)]), 1]]

        expansion_term=[[np.array([0 for coupling in range(self.couplings_amount)]), 1]]

        for order in range(1, self.max_PT_order+1):

            expansion_term=[[C[0], C[1]*(alpha-order+1)/order] for C in self.PT_cutoff(C_series_mult(expansion_term,X))]
            N=C_series_add(N,expansion_term)

            self.log(f'added to N: {expansion_term}\n')
            self.log(f'new N: {N}\n')

        N=sorted(N, key= lambda C: np.sum(C[0]))

        return( self.PT_cutoff( N) )


    def normalize_C_dictionary(self, the_C_dictionary):
    
    #     to-be class method:
    # test by running it directly or as a part of a larger procedure

        the_C_list=[[np.array(eval(key)), the_C_dictionary[key]] for key in the_C_dictionary]

        Z=self.PT_cutoff(self.C_series_to_Z(the_C_list))

        N=self.N_from_Z(Z)

        normalized_C_list=self.PT_cutoff(C_series_mult(the_C_list, N))

        normalized_C_dictionary={str(a_C[0].tolist()): a_C[1] for a_C in normalized_C_list}



        return(normalized_C_dictionary)

    def eq_adapt_C_series(self, C_series):
    
    #     to-be class method:
    # test by running it directly or as a part of a larger procedure

        '''
        C_series is assumed to have raw form ({'[0,0,0]': 1, ...}), 
        but the output is in the equation-ready form ({'K = {[K]}, s = {s}': ...})
        '''

        adapted_C_series = dict()

        for K in C_series:

            if K == str([0 for coupling in range(self.couplings_amount)]):

                continue

            if self.PT_mode=='scalar':

                C_K_label=f'K = {[sum(eval(K))]}, s = {self.s_of_k(eval(K))}'

            elif self.PT_mode=='vector':

                C_K_label=f'K = {K}, s = {self.s_of_k(eval(K))}'

            else:

                raise Exception('Unknown PT_mode!')

            if C_K_label in adapted_C_series:

                adapted_C_series[C_K_label] += C_series[K]

            else:

                adapted_C_series.update({C_K_label: C_series[K]})

            self.log(f'included the value of {C_K_label} type into the series: {C_series[K]}')

        return (adapted_C_series)
    
    
    ### Wavefunction representation

    def wavefunction_from_PT_series (self, normalized_C_dictionary, J):
    
    # to-be class method:
    # test by running it directly or as a part of a larger procedure

        computational_states=[ [[np.array([1 if s==label else 0 for label in range(2)]) for s in s_string], coef] for s_string, coef in [[self.s_of_k(eval(K)), normalized_C_dictionary[K]*J**sum(eval(K))] for K in normalized_C_dictionary]]

        wavefunction=sum([coef*reduce(np.kron, computational_state) for computational_state, coef in self.computational_states])

        wavefunction=wavefunction/np.linalg.norm(wavefunction)

        return(wavefunction)

        
    def theta_to_k_correspondence(self, theta_k):
    #     to-be class method:
    # move it to the class body
    # test by running it directly or as a part of a larger procedure

        return(k_to_generator(theta_k[1])==self.theta_to_generator(theta_k[0]))

    def f_theta_PT_filter(self, f_theta):
    #     to-be class method:
    # move it to the class body
    # test by running it directly or as a part of a larger procedure    

        return(sum([sum(theta_k[1]) for theta_k in f_theta])<=self.max_PT_order)


    def theta_k_filter(self, theta_k):

    #     to-be class method:
    # move it to the class body
    # test by running it directly or as a part of a larger procedure        

        if self.PT_mode=='vector':
            return(odd_count_equals_one(theta_k) and self.theta_to_k_correspondence(theta_k))

        elif self.PT_mode=='scalar':
            return(odd_count_equals_one(theta_k))

        else:
            raise Exception('Unknown PT_mode!')


    
    ## f-analysis

    ### $\vec{K}(f)$, $\vec{N}(f)$ and $\Theta(f)$

    def f_theta_to_pauli_strings(self, f_theta):
    # to-be class method:
    # test by running it directly or as a part of a larger procedure    


        return(self.pauli_strings_from_generator_list([self.theta_to_generator(k_theta[0]) for k_theta in f_theta]))

    def comp_state_from_f_theta(self, f_theta):
    # to-be class method:
    # test by running it directly or as a part of a larger procedure

        generators_action=threaded_pauli_strings_action(self.f_theta_to_pauli_strings(f_theta),
                                                        [[0 for spin_number in range(self.spin_amount)], 1j**len(f_theta)])

        return(generators_action)

    def product_function(self, theta_ks, f_theta, theta_k_variable_list):

    #     to-be class method:
    # move it to the class body
    # test by running it directly or as a part of a larger procedure

        theta_k_indices=[theta_k_variable_list.index(k_theta) for k_theta in f_theta]

        prefactor=self.comp_state_from_f_theta(f_theta)[1]/multiplicities_factorials(f_theta)

        return(reduce(mul, [theta_ks[index] for index in theta_k_indices] )* prefactor)

    def theta_function(self, theta_ks, fs_theta, theta_k_variable_list):

#     to-be class method:
# move it to the class body
# test by running it directly or as a part of a larger procedure
    
        return(sum(self.product_function(theta_ks, f, theta_k_variable_list) for f in fs_theta))
# Section 5 of the jupyter notebook


    ### f list generation


    def f_theta_set_function(self):

    #     to-be class method:
    # move it to the class body
    # test by running it directly or as a part of a larger procedure

        stopwatch=time.time()

        PT_orders_for_theta=self.list_of_Ks_from_PT_order()

        self.log(f'PT_orders_for_theta: {PT_orders_for_theta}')

        theta_k_set=set()

        for theta in range(self.number_of_thetas):
            theta_k_set.update([(theta, tuple(a_K)) for a_K in PT_orders_for_theta ])

        self.log(f'theta_k_set: {theta_k_set}')

    #     For TUCC applied to TFIM, it makes sense to remove all even PT orders from theta expansion
    #     theta_k_set=set(filter(lambda theta_k: self.theta_k_filter(theta_k, PT_mode), theta_k_set))

        self.log(f'filtered theta_k_set: {theta_k_set}')

        all_f_thetas=set()

        new_power_f_thetas=set((theta_k,) for theta_k in theta_k_set)

        self.log(f'new_power_f_thetas: {new_power_f_thetas}')

        all_f_thetas.update(new_power_f_thetas)

        for theta_power in range(1, self.max_PT_order+1):   

            '''
            Simply listing all potential next-power terms
            '''

            new_power_f_thetas=set(a_product[0]+(a_product[1],) for a_product 
                                   in set(itertools.product(new_power_f_thetas,theta_k_set)) )

            self.log(f'new_power_f_thetas: {new_power_f_thetas}')


            '''
            Sorting the f_thetas: to avoid dublicates and to order terms for the T-action
            '''

            new_power_f_thetas=set(tuple(sorted(list(f_theta), key=lambda theta_k: theta_k[0])) for f_theta in new_power_f_thetas)

            self.log( f'sorted new_power_f_thetas: {new_power_f_thetas}')

            '''
            Filtering out the higher order terms
            '''

            new_power_f_thetas=set(filter(lambda f_theta: self.f_theta_PT_filter( f_theta), new_power_f_thetas))

            self.log(f'filtered new_power_f_thetas: {new_power_f_thetas}')

            all_f_thetas.update(new_power_f_thetas)

            self.log(f'all_f_thetas: {all_f_thetas}')

        self.log(f'f_theta_set evaluation completed, time elapsed: {time.time()-stopwatch}')

        return(list(theta_k_set), all_f_thetas)



    def f_theta_to_K_s_dict(self, the_f_theta_set):

    #     to-be class method:
    # move it to the class body
    # test by running it directly or as a part of a larger procedure

        f_theta_K_s_dict=dict()

        for f_theta in the_f_theta_set:

            K_for_f_theta=PT_order_from_f_theta(f_theta)

            s_for_f_theta=tuple(self.comp_state_from_f_theta(f_theta)[0])

            label_for_f_theta=f'K = {str(PT_order_from_f_theta(f_theta))}, s = {self.comp_state_from_f_theta(f_theta)[0]}'

            if label_for_f_theta in f_theta_K_s_dict:

                f_theta_K_s_dict[label_for_f_theta].update({f_theta})

            else:

                f_theta_K_s_dict.update({label_for_f_theta : {f_theta}})

        for K_s in f_theta_K_s_dict:
            self.log( f'{K_s}: {f_theta_K_s_dict[K_s]}')    

        return(f_theta_K_s_dict)

    
    def equation_initialize (self, equation_size_statistics=True):
        
        stopwatch=time.time()
        
        C_series=self.normalize_C_dictionary(self.unnormalized_C_dictionary())

        eq_adapted_C_series=self.eq_adapt_C_series(C_series)

        theta_k_variable_list, f_theta_set_for_eq = self.f_theta_set_function()

        f_dict_for_eq = self.f_theta_to_K_s_dict(f_theta_set_for_eq)

        list_of_equations=list(eq_adapted_C_series)

        for key in f_dict_for_eq:

            if key in eq_adapted_C_series:

                pass

            else:

                eq_adapted_C_series[key]=0

        for key in eq_adapted_C_series:

            if key in f_dict_for_eq:

                pass

            else:

                f_dict_for_eq[key]={}      
        

        self.eq_adapted_C_series=eq_adapted_C_series

        self.theta_k_variable_list=theta_k_variable_list

        self.f_theta_set_for_eq=f_theta_set_for_eq

        self.f_dict_for_eq=f_dict_for_eq

        self.list_of_equations=list_of_equations
        
        self.log(f'Equations initialization complete in {time.time()-stopwatch} seconds', True)
        
        self.log(f'The number of variables is: {len(self.theta_k_variable_list)}', equation_size_statistics)
        
        
        
        no_of_equations=len(self.list_of_equations)
        
        avg_num_of_terms=sum([len(self.f_dict_for_eq[key]) for key in self.f_dict_for_eq])/no_of_equations
                                                                                               
        avg_nonlinearity=sum([len(term) for key in self.f_dict_for_eq for term in self.f_dict_for_eq[key]])/(avg_num_of_terms*no_of_equations)                   
                                                                                               
        self.log(f'The number of equations is: {no_of_equations}', equation_size_statistics)
        
        self.log(f'The avg number of terms per equation is: {avg_num_of_terms}', equation_size_statistics)
        
        self.log(f'The avg nonlinearity is: {avg_nonlinearity}', equation_size_statistics)

        
                                                                                              
        return(None)
    
        
    def equation_system(self, theta_k_values, option='non-weighted'):

#     For some reason ignoring the J**K weight during the optimization is sometimes beneficial from the perspective of the
#     _weighted_ function. I will comment this part out for now, but we may look at it in closer detail later.



        if option=='non-weighted':
    
            list_of_eq_values=[np.abs(self.theta_function(theta_k_values, self.f_dict_for_eq[equation], self.theta_k_variable_list)
                         - self.eq_adapted_C_series[equation])  
                               for equation in self.list_of_equations]
        
        elif option=='weighted':
            
            list_of_eq_values=[np.abs(self.theta_function(theta_k_values, self.f_dict_for_eq[equation], self.theta_k_variable_list)
                     - self.eq_adapted_C_series[equation])  * (self.J**eval(re.search('K = \[(.+?)\]',equation).group(1)) )
                           for equation in self.list_of_equations]
        else:
            
            raise Exception('Option not recognised!')

    #   "eval(re.search('K = \[(.+?)\]',equation).group(1))" is the value of K in the equation

        

        if len(self.list_of_equations)<len(self.theta_k_variable_list):
            list_of_eq_values+=[0 for additional_equation in range(len(self.theta_k_variable_list)-len(self.list_of_equations))]

            self.log(f'list_of_eq_values={list_of_eq_values}')

        return(np.array(list_of_eq_values, dtype='float'))
    
    def theta_k_to_theta_values(self):
        
        
        #assumes that self.theta_k_values already exists
        
        theta_values=[0 for theta_label in range(len(self.unitary_generators))]
        
        for contribution_number in range(len(self.theta_k_variable_list)):

            contribution_label=self.theta_k_variable_list[contribution_number]

            theta_index=contribution_label[0]

            degree_of_PT=contribution_label[1][0]

            contribution_value=self.theta_k_values[contribution_number]

            theta_values[theta_index]+=contribution_value*self.J**degree_of_PT
        
        return(theta_values)
        
    
    def equation_solving (self, option='non-weighted'):
        
        stopwatch=time.time()

        self.theta_k_values=optimize.minimize(lambda theta_ks: sum(self.equation_system(theta_ks, option)**2), np.array([0 for k_theta in self.theta_k_variable_list]), method='SLSQP', bounds=tuple((-1/self.J, 1/self.J) for theta_k in range( len(self.theta_k_variable_list) )) ).x
        
#         self.theta_k_values=optimize.leastsq(self.equation_system, np.array([0 for k_theta in self.theta_k_variable_list]))[0]
        
        
        self.log(f'The bounds are: {tuple((-1/self.J, 1/self.J) for theta_k in range( len(self.theta_k_variable_list) ))}', False)

        self.log(f'Equations solved in: {time.time()-stopwatch} seconds,\n The equations are solved with precision:  {sum(self.equation_system(self.theta_k_values, option)**2)} \n, The solutions are returned and also stored in self.theta_values', True)
        
        self.theta_values=self.theta_k_to_theta_values()
        
        self.precision=sum(self.equation_system(self.theta_k_values,option)**2)
        
        return(self.theta_values, self.precision)
    
# Tools and functions

## General tools

### Combinatorics 

def partitioning (n):
    
    #confirmed external function
    
    '''
    Creates all possible partitions of an integer n, as a list of lists
    '''
    partitions=[[n]]
    

    
    while partitions[-1][0]>1:
    
        partition=partitions[-1]

        '''
        Finding the rightmost non-one, reducing it by one:
        '''
            
        for k in range(len(partition)):

            if partition[k]>1:
                hit_k=k

            else:
                break
        
        hit_p=partition[hit_k]-1

        rest=sum(partition[hit_k+1:])+1
        
        '''
        Stacking up the rest:
        '''
        
        assemble=[hit_p for i in range(rest//hit_p) ]

        if rest-sum(assemble)!=0:
            assemble+=[rest-sum(assemble)]
        
        '''
        Gathering up the new partition:
        '''
        
        partitions+=[partition[:hit_k]+[hit_p]+assemble]
        
    

    return(partitions)

def fill_the_list(the_list, full_length, filler=0):
    
    #confirmed external function
    
    if len(the_list)<=full_length:
        return(the_list+[0 for iterator in range(full_length-len(the_list))])
    else:
        return(None)

### k vectors, couplings and generators


def k_to_couplings(k):
#     !!!!!! TFIM-only thing
#     confirmed non-class method
    
    '''
    Assumes that the couplings are commuting Paulis, and implicitly that only one coupling i has k_i%2==1
    
    Returns a list of coupling labels as integers from 1 to N_c
    '''
    
    k=list(k)
    
    return([coupling_label for coupling_label in range(len(k)) if k[coupling_label]%2==1])



def k_to_generator(k):
    
#     confirmed non-class method
#     seems to be for a specific coupling type!

    
    '''
    Assuming odd-even logic in the couplings. 
    k is either list or a numpy array
    '''
    
    k=list(k)
    
    for element in k:
        if element%2==1:
            return k.index(element)
        
    print('Didn"t find a generator for k!')
    
    return(None)




### Pauli action

def operator_from_pauli_string(string):
    
#     confirmed non-class method
    
    single_qubit_operator_list=[]
    
    for element in string:
        if element=='X':
            single_qubit_operator_list+=[np.array([[0,1],[1,0]])]
        elif element=='Y':
            single_qubit_operator_list+=[np.array([[0,-1j],[0,1j]])]
        elif element=='Z':
            single_qubit_operator_list+=[np.array([[1,0],[0,-1]])]
        elif element=='I':
            single_qubit_operator_list+=[np.array([[1,0],[0, 1]])] 
        else:
            raise('Non-pauli input')
    
    return(reduce(np.kron, single_qubit_operator_list))

def single_pauli_action(pauli, spin):
    
#     confirmed non-class method
    
    if pauli=='X':
        return(np.mod(spin+1,2), 1)
    elif pauli=='Y':
        return(np.mod(spin+1,2), 1j*(-1)**spin)
    elif pauli=='Z':
        return(spin, (-1)**spin)
    elif pauli=='I':
        return(spin, 1)
    else:
        print('wrong pauli!')
        return(None)

def pauli_string_action(pauli_string, spins_and_prefactor):
    
    '''
    Given a pauli_string (label) and a computation basis state+prefactor, 
    returns a new computational state with a new prefactor
    
    spins_and_prefactors=[integer list, complex number]
    '''
#     confirmed non-class method
    
    spins=spins_and_prefactor[0]
    
    assert len(pauli_string)==len(spins)
    
    new_spins_and_prefactor=[single_pauli_action(pauli_string[spin_number], spins[spin_number]) for spin_number in range(len(spins))]
    
    return([[element[0] for element in new_spins_and_prefactor], 
                         spins_and_prefactor[1]*reduce(lambda x,y: x*y, [element[1] for element in new_spins_and_prefactor])])


def threaded_pauli_strings_action(pauli_strings,spins_and_prefactor):
    
#     confirmed non-class method
    
    return(reduce(lambda s_and_p, p_str: pauli_string_action(p_str, s_and_p), [spins_and_prefactor]+pauli_strings) )

## PT series

### Dyson elementaries

def E0_of_s(s):
    
#     confirmed non-class method

    return( sum([-1*(-1)**spin_value for spin_value in s]))



### C-series manipulations

def C_mult(C1, C2):

#     confirmed non-class method
    
    return([C1[0]+C2[0], C1[1]*C2[1]])



def C_series_add(C_series_1, C_series_2):
    
#     confirmed non-class method
    
    the_sum=C_series_1+C_series_2
    
    K_set={tuple(C[0]) for C in the_sum}

    
#     log(f'add={add}')
#     log(f'K_set={K_set}')
    
    the_sum=[[np.array(K), sum([C[1] for C in the_sum if tuple(C[0])==K])] for K in K_set]
    
    the_sum=sorted(the_sum, key=lambda C: np.sum(C[0]))
    
    return(the_sum)

def C_series_mult(C_series_1, C_series_2):
    
#     confirmed non-class method
    
    the_product=[[C1[0]+C2[0], C1[1]*C2[1]] for C1, C2 in list(itertools.product(C_series_1, C_series_2))]
    
    K_set={tuple(C[0]) for C in the_product}

    
#     log(f'mult={mult}')
#     log(f'K_set={K_set}')
    
    the_product=[[np.array(K), sum([C[1] for C in the_product if tuple(C[0])==K])] for K in K_set]
    
    the_product=sorted(the_product, key=lambda C: np.sum(C[0]))
    
    return(the_product)





## f-analysis

### $\vec{K}(f)$, $\vec{N}(f)$ and $\Theta(f)$

def PT_order_from_f_theta(f_theta):
    
#     confirmed non-class method
    
    return(list(sum([np.array(k_theta[1]) for k_theta in f_theta])))



def multiplicities_factorials(f_theta):

    #     confirmed non-class method

    counts=[]
    for element in f_theta:
        counts+=[f_theta.count(element)]
        f_theta=list(filter(lambda element_prime: element_prime!=element, f_theta))

    return(reduce(mul, [math.factorial(count) for count in counts]) )
    


### Combinatoric tools

def odd_count(the_list):
    
    #     confirmed non-class method
    
    return(len([element for element in the_list if element%2==1]))

def odd_count_equals_one(theta_k):
    
    #     confirmed non-class method
    
    return(odd_count(theta_k[1])==1)

### f list generation






def f_theta_to_K_dict(the_f_theta_set):
    
    #confirmed non-class method
    
    f_theta_K_dict=dict()
    
    for f_theta in the_f_theta_set:
        
        label_for_f_theta = f'K = {str(PT_order_from_f_theta(f_theta))}'
        
        if label_for_f_theta in f_theta_K_dict:
            
            f_theta_K_dict[label_for_f_theta].update({f_theta})
            
        else:
            
            f_theta_K_dict.update({label_for_f_theta : {f_theta}})
       
        
    return(f_theta_K_dict)


//
//  censoredsplitmerge
//
//   Simulate the coupling used by Schramm in the paper Compositions of random transpositions [ Israel J. Math. 147 221–243],
//  but generalised so that a proposed merge of two blocks is only accepted after flipping a biased coin (probability beta_m in the
//  code below), as in Mayer-Wolf, Zeitouni, Zerner.
//
//  This modified process has Poisson-Dirichlet(1/beta_m) invariant distribution.
//
//  I didn't implement rejection of splits, as Schramm's proof works in this case.
//
//  PW


// Dumps output to /tmp/censoredsplitmerge.csv
// I load the data into R like this:
// S<-read.csv("/tmp/censoredsplitmerge.csv", header=F)

#include <iostream>
#include <numeric>
#include <iterator>
#include <algorithm>
#include <deque>
#include <random>
#include <fstream>
#include <cassert>

// merge probability
const double beta_m = 0.6;

// number of updates to perform
const int n_steps = 250000;

// unmatched parts in Y, Z, and YZ = common matched parts
std::deque<double> Y, Z, YZ;

// keep track of matched mass (this could be computed dynamically but it is more efficient just to update it on the fly)
double YZ_mass;

// random number generator (Mersenne Twister with default seed)
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> Unif(0, 1);

// for the merge acceptance
std::bernoulli_distribution accept_merge_coin(beta_m);


// takes a deque of doubles, looks at the cumulative sum and returns
// an iterator pointing at the last element not exceeding u
// If u exceeds the sum of all elements then partition.end() is returned
// partition should be non-empty
std::deque<double>::iterator find_part(double u, std::deque<double> & partition) {
    double cum_sum = 0;
    std::deque<double>::iterator part;
    for (part = partition.begin(); part != partition.end(); ++part ) {
        cum_sum += *part;
        if (cum_sum > u) break;
    }
    return part;
}


// ** The followign functions are only for a given deque;
//   i.e. they don't move mass to/from the matched deques ** //

// merges the first part into a given part
void merge_part(std::deque<double> & partition, std::deque<double>::iterator part) {
    *part += *partition.begin();
    partition.pop_front();
}

// splits the part at the beginning into a part of size v and whatever is left.
// assert() the size of the first part is larger than v.
void split_first(std::deque<double> & partition, double v) {
    assert(*partition.begin() >= v);
    partition.push_back(*partition.begin() - v);
    *partition.begin() = v;
}

// find the biggest part; it is ok for the partitions to be empty here
double biggest_part(std::deque<double> &partition){
    std::deque<double>::iterator maxpart = std::max_element(partition.begin(),partition.end());
    if (maxpart == partition.end()) {
        return 0.0;
    } else return *std::max_element(partition.begin(),partition.end());
}

double second_biggest_part(std::deque<double> &partition){
    std::deque<double>::iterator maxpart = std::max_element(partition.begin(),partition.end());
    if (maxpart == partition.end()) {
        return 0.0;
    } else {
        double tmp = *maxpart;
        *maxpart = 0;
        double secondpart = *(std::max_element(partition.begin(),partition.end()));
        *maxpart = tmp;
        return secondpart;
    }
}

int main(int argc, const char * argv[])
{
    //  ******** initial conditions *********
    Y.push_back(1.0);
    Z.push_back(1.0);
    //Z.push_back(0.5);

    
    // for initial distribution just run the split merge processes independently for a while
    
    for (int j = 0; j <5000 ; j++) {
        // generate U, V = Unif(0,1)
        double U = Unif(gen), V = Unif(gen);
        
        // move chosen part to the front
        iter_swap(Y.begin(), find_part(U,Y));
        
        // find the second part
        std::deque<double>::iterator pVY = find_part(V,Y);
        
        if (pVY == Y.begin()) { // we chose same part twice; split it at V
            split_first(Y, V);
        } else if(accept_merge_coin(gen)) { // else it is a merge proposal so flip the coin and merge if necessary
            merge_part(Y, pVY);
        }
    }
    
    // do same for Z
    for (int j = 0; j <5000 ; j++) {
        // generate U, V = Unif(0,1)
        double U = Unif(gen), V = Unif(gen);
        
        // move chosen part to the front
        iter_swap(Z.begin(), find_part(U,Z));
        
        // find the second part
        std::deque<double>::iterator pVZ = find_part(V,Z);
        
        if (pVZ == Z.begin()) { // we chose same part twice; split it at V
            split_first(Z, V);
        } else if(accept_merge_coin(gen)) { // else it is a merge proposal so flip the coin and merge if necessary
            merge_part(Z, pVZ);
        }
    }
    
    
    // ******* coupled updates **********
    
    // initially no matched mass (note can setup initially matched mass above)
    // if (YZ.empty()) YZ.push_back(0.0);

    // variable for the matched mass (i.e. sum of all elements in YZ
    YZ_mass = std::accumulate(YZ.begin(), YZ.end(),0.0);
    
    
    // log file for results in CSV format (for importing into R for analysis)
    std::ofstream log_file;
    log_file.open("/tmp/censoredsplitmerge.csv");
    
    for (int j = 0; j < n_steps ; j++)
    {
        
        // matched mass, largest matched, Y_1, Y_2, Z_1, Z_2,
        log_file << YZ_mass << ", " <<  biggest_part(YZ) << ","  << biggest_part(Y) << ", " << second_biggest_part(Y) << "," <<  biggest_part(Z) << "," << second_biggest_part(Z) << "," << Y.size() << "," << Z.size() << "," << YZ.size() << ",";
    
        bool generated_matched_mass = false;
        
        bool prop_SM = false;
        
        // sort the deques so blocks of similar size are aligned better
        std::sort (Y.begin(), Y.end());
        std::sort (Z.begin(), Z.end());
        
        // generate U, V = Unif(0,1)
        double U = Unif(gen), V = Unif(gen);
        
        if (U < YZ_mass) {  // **** first part chosen from matched mass
            
            assert(!YZ.empty()); // [sanity check]
            
            // find the corresponding part in the matched mass and move it to the front
            iter_swap(YZ.begin(), find_part(U,YZ));
            
            // [The chosen part is now at YZ.begin()]
            
            if (V < YZ_mass) { // **** both parts from matched mass
                std::deque<double>::iterator p = find_part(V,YZ);
                
                if (p == YZ.begin()) { // split at V
                    split_first(YZ, V);
                }
                else if(accept_merge_coin(gen)){ // flip the biased coin, and merge if heads
                    merge_part(YZ,p);
                }
            } else { // **** first part matched, second unmatched
                // only possibility is to merge in both Y,Z and lose the matched mass
                if(accept_merge_coin(gen)) {
                    assert(!Y.empty() && !Z.empty());
                    
                    std::deque<double>::iterator pY = find_part(V - YZ_mass,Y), pZ = find_part(V - YZ_mass,Z);
                    
                    // merge the parts
                    *pY += *YZ.begin();
                    *pZ += *YZ.begin();
                    
                    // bye bye matched mass :(
                    YZ_mass -= *YZ.begin();
                    YZ.pop_front();
                    
                    // keep YZ non-empty
                    if (YZ.empty()) YZ.push_back(0.0);
                }
            } // (second choice using V)
        }
        else {  // ****** first part chosen from unmatched mass
            assert(!Y.empty() && !Z.empty());
            
            // move the chosen parts to the beginning of the deque
            iter_swap(Y.begin(), find_part(U - YZ_mass,Y));
            iter_swap(Z.begin(), find_part(U - YZ_mass,Z));
            
            
            if (V < YZ_mass) { // ***** second part chosen from matched mass
                assert(!YZ.empty());
                
                // only possibility is a merge in both Y and Z with a matched part
                if(accept_merge_coin(gen)) {
                    std::deque<double>::iterator p = find_part(V,YZ);
                    *Y.begin() += *p;
                    *Z.begin() += *p;
                    // remove the matched part
                    YZ_mass -= *p;
                    YZ.erase(p);
                }
            }
            else { // ***** both parts chosen from unmatched mass;
                // situation is more complicated as we may have a split in Y and a merge in Z (or vice versa)
                
                // Let us find the corresponding parts
                std::deque<double>::iterator pY = find_part(V - YZ_mass,Y), pZ = find_part(V - YZ_mass,Z);
                
                // We sample the coin in advance, so that it can be shared between Y and Z if merges are proposed
                bool accept_merge = accept_merge_coin(gen);
                
                if (pY == Y.begin()) { // split first block
                    assert(*pY > V - YZ_mass);
                    split_first(Y,V - YZ_mass);
                    prop_SM = true;
                }
                else {
                    prop_SM = false;
                    if (accept_merge) { // merge
                        merge_part(Y, pY);
                    }
                }
                
                if (pZ == Z.begin()) { // split first block
                    assert(*pZ > V - YZ_mass);
                    split_first(Z,V - YZ_mass);
                    prop_SM = !prop_SM;
                }
                else {
                    if (accept_merge) {
                        merge_part(Z, pZ);
                    }
                }
                
                if (*Y.begin() == *Z.begin()) {
                    // generate new matched mass
                    YZ.push_back(*Y.begin());
                    YZ_mass += *Y.begin();
                    
                    // remove from unmatched
                    Y.pop_front();
                    Z.pop_front();
                    
                    generated_matched_mass = true;
                }
                
            }
        }
        
                
        if (prop_SM) {
            log_file << "1";
        } else {
            log_file << "0";
        }

        log_file << std::endl;
        
    } // end update step
    
    log_file.close();
    

    return 0;
}

//    double Y_mass, Z_mass, YZ_mass_check;
//
//    Y_mass = std::accumulate(Y.begin(), Y.end(),0.0);
//    Z_mass = std::accumulate(Z.begin(), Z.end(),0.0);
//    YZ_mass_check = std::accumulate(YZ.begin(), YZ.end(),0.0);
//    //
//    //        std::cout << "Y = {";
//    //        std::copy(Y.begin(), Y.end(), std::ostream_iterator<double>(std::cout, ", "));
//    //        std::cout << "}" << ", mass = " << Y_mass << std::endl;
//    //
//    //
//    //        std::cout << "Z = {";
//    //        std::copy(Z.begin(), Z.end(), std::ostream_iterator<double>(std::cout, ", "));
//    //        std::cout << "}" << ", mass = " << Z_mass << std::endl;
//    //
//    //        std::cout << "YZ = {";
//    //        std::copy(YZ.begin(), YZ.end(), std::ostream_iterator<double>(std::cout, ", "));
//    //        std::cout << "}" << ", mass = " << YZ_mass_check << std::endl;
//    //
//    //
//    std::cout << YZ_mass_check - YZ_mass << ", " << Y_mass - Z_mass << ", " << 1 - Z_mass - YZ_mass << std::endl;



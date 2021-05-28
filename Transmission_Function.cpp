/*
 * This script contains the Rcpp code for the nested "for" loop that controls
 * interactions between agents.
 */



#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::interfaces(r, cpp)]]

// Create a constant integer value for integer arithmetic
// int CONSTANT = 1000;

// Here is the function that propagates infection based on agent interactions.

// [[Rcpp::export]]
Rcpp::List propogate(Rcpp::List agents, Rcpp::IntegerMatrix contact_matrix, Rcpp::NumericVector initial_set)
{

  // Initialize necessary variables
  int len = 0, num_hcws = 0, num_pats = 0, val = 0, send_disease_stat = 0;
  int receive_disease_stat = 0, num_contacts = 0, new_infects = 0;
  int hcw_to_pat = 0, pat_to_hcw = 0, secondary_hcw_infected = 0, id = 0;
  int secondary_patient_infected = 0, Rnot = 0, val1 = 0, val2 = 0, pat_to_pat = 0, hcw_to_hcw = 0;
  double infect = 0.0, sus = 0.0, random_double_one = 0.0, random_double_two = 0.0;

  // Create lists to hold agents
  Rcpp::List ids, infects, suscs, types, disease_stats;

  // Rcpp::Rcout << "Rcpp is running." << std::endl;

  // Extract lists of different agent characteristics
  ids           = agents("id");
  types         = agents("type");
  infects       = agents("infect");
  suscs         = agents("sus");
  disease_stats = agents("disease_stat");

  len = ids.length();
  // Rcpp::Rcout << "Number of agents is " << len << std::endl;

  // For loop sees how many agents of each type there are
  for (int i = 0; i < len; i++)
  {

    val = types[i];

    if (val == 1)
    {
      num_pats++;
    }

    if (val == 2)
    {
      num_hcws++;
    }


  }

  // Rcpp::Rcout << "Number of HCWs is " << num_hcws << std::endl;
  // Rcpp::Rcout << "Number of patients is " << num_pats << std::endl;

  // For loop goes through all the infecting agents
  for (int i = 0; i < len; i++)
  {

    // Find infectivity of sending agent
    infect = infects[i];

    // Find disease status of sending agent
    send_disease_stat = disease_stats[i];

    if (infect == 0 || send_disease_stat != 1)
    {
      continue;
    }

    // For loop goes through all the susceptible agents
    for (int j = 0; j < len; j++)
    {

      // Find susceptibility, disease status, and number of contacts  of receiving agent
      sus                  = suscs[j];
      receive_disease_stat = disease_stats[j];
      num_contacts         = contact_matrix(i, j);

      if (sus == 0 || receive_disease_stat != 0 || num_contacts == 0)
      {
        continue;
      }
      else
      {

       // Check if contacts are less than one for timestemp
       if (num_contacts < 1)
       {
         if (R::runif(0, 1) < num_contacts)
         {
            num_contacts = 1;
         }
         else
         {
           continue;
         }
       }

        // For loop draws a random number for each contact to see whether susceptible agent gets infected.
        for (int k = 0; k < num_contacts; k++)
        {

          random_double_one = R::runif(0, 1);
          random_double_two = R::runif(0, 1);
          if (random_double_one <= infect && random_double_two <= sus)
          {
            // If infected, disease status becomes 2
            disease_stats[j] = 2;
            new_infects++;

            // Update R0 if necessary
            id = ids[i];

            if (std::count(initial_set.begin(), initial_set.end(), id))
            {
              Rnot++;
            }

            // Check who is infecting whom.
            val1 = types[i];
            val2 = types[j];

            if (val2 == 1)
            {
              secondary_patient_infected++;
            }
            else if (val2 == 2)
            {
              secondary_hcw_infected++;
            }

            if (val1 == 1 && val2 == 1)
            {
              pat_to_pat++;
            }
            else if(val1 == 2 && val2 == 1)
            {
              hcw_to_pat++;
            }
            else if (val2 == 2 && val1 == 2)
            {
              hcw_to_hcw++;
            }
            else if (val1 == 1 && val2 == 2)
            {
              pat_to_hcw++;
            }

            break;
          }

        }
      }

    }
  }

  disease_stats.push_back(new_infects);
  disease_stats.push_back(hcw_to_pat);
  disease_stats.push_back(pat_to_hcw);
  disease_stats.push_back(pat_to_pat);
  disease_stats.push_back(hcw_to_hcw);
  disease_stats.push_back(secondary_patient_infected);
  disease_stats.push_back(secondary_hcw_infected);
  disease_stats.push_back(Rnot);

  return disease_stats;

}

/*
###############################################################################
# HCC KHK Metabolism Model
###############################################################################
*/

#include "./custom.h"

void create_cell_types( void )
{
    SeedRandom( parameters.ints("random_seed") ); 
    
    initialize_default_cell_definition();
    cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
    
    cell_defaults.functions.update_phenotype = NULL; 
    cell_defaults.functions.custom_cell_rule = NULL; 
    cell_defaults.functions.contact_function = NULL; 
    
    cell_defaults.name = "HCC cell"; 
    cell_defaults.type = 0; 
    
    // Add custom data
    cell_defaults.custom_data.add_variable( "KHK_expression" , "dimensionless", 17.0 );
    cell_defaults.custom_data.add_variable( "glucose_uptake_rate" , "mg/min", 0.0 );
    cell_defaults.custom_data.add_variable( "fructose_uptake_rate" , "mg/min", 0.0 );
    cell_defaults.custom_data.add_variable( "ATP_level" , "dimensionless", 1.0 );
    
    // Phenotype
    cell_defaults.phenotype.cycle.data.transition_rate(0,0) = 0.04/24.0;
    cell_defaults.functions.update_phenotype = tumor_cell_phenotype; 
    
    build_cell_definitions_maps(); 
    display_cell_definitions( std::cout ); 
    
    return; 
}

void setup_microenvironment( void )
{
    default_microenvironment_options.X_range = {-500, 500}; 
    default_microenvironment_options.Y_range = {-500, 500}; 
    default_microenvironment_options.Z_range = {-50, 50}; 
    
    default_microenvironment_options.simulate_2D = false; 
    
    microenvironment.add_density( "glucose", "mM" , 1e5, 0.01 ); 
    microenvironment.add_density( "fructose", "mM" , 1e5, 0.001 ); 
    
    default_microenvironment_options.outer_Dirichlet_conditions = true;
    default_microenvironment_options.Dirichlet_condition_vector[0] = 5.5;
    default_microenvironment_options.Dirichlet_condition_vector[1] = 0.5;
    
    initialize_microenvironment(); 
    
    return; 
}

void setup_tissue( void )
{
    Cell* pC = create_cell(); 
    pC->assign_position( 0.0, 0.0, 0.0 );
    
    double khk_value = parameters.doubles("KHK_expression");
    pC->custom_data["KHK_expression"] = khk_value;
    
    std::cout << "Initialized tumor with KHK = " << khk_value << std::endl;
    
    return; 
}

void tumor_cell_phenotype( Cell* pCell, Phenotype& phenotype, double dt )
{
    static int glucose_index = microenvironment.find_density_index( "glucose" );
    static int fructose_index = microenvironment.find_density_index( "fructose" );
    
    double KHK = pCell->custom_data["KHK_expression"];
    double glucose = pCell->nearest_density_vector()[glucose_index];
    double fructose = pCell->nearest_density_vector()[fructose_index];
    
    double Km_glucose = 7.0;
    double Vmax_glucose = 0.02;
    double glucose_uptake = (Vmax_glucose * glucose) / (Km_glucose + glucose);
    
    double Km_fructose = 6.0;
    double Vmax_fructose_base = 0.015;
    double KHK_normalized = std::max(0.0, std::min(1.0, (KHK - 14.0) / 5.0));
    double Vmax_fructose = Vmax_fructose_base * (0.5 + 1.5 * KHK_normalized);
    double fructose_uptake = (Vmax_fructose * fructose) / (Km_fructose + fructose);
    
    pCell->custom_data["glucose_uptake_rate"] = glucose_uptake;
    pCell->custom_data["fructose_uptake_rate"] = fructose_uptake;
    
    double ATP = (glucose_uptake * 30.0 + fructose_uptake * 25.0) / 50.0;
    pCell->custom_data["ATP_level"] = ATP;
    
    double ATP_factor = std::max(0.0, std::min(1.0, (ATP - 0.5) / 0.5));
    double substrate_factor = glucose / (glucose + 2.0);
    phenotype.cycle.data.transition_rate(0,0) = (0.04/24.0) * ATP_factor * substrate_factor;
    
    phenotype.secretion.uptake_rates[glucose_index] = glucose_uptake * 60.0;
    phenotype.secretion.uptake_rates[fructose_index] = fructose_uptake * 60.0;
    
    if (ATP < 0.2) { pCell->start_death(0); }
    
    return;
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{
    std::vector<std::string> output = {"black", "black", "black", "black"};
    double KHK = pCell->custom_data["KHK_expression"];
    
    if (KHK < 16.0) {
        output[0] = "blue"; 
        output[2] = "darkblue";
    } else if (KHK < 18.0) {
        output[0] = "yellow"; 
        output[2] = "orange";
    } else {
        output[0] = "red"; 
        output[2] = "darkred";
    }
    
    return output;
}

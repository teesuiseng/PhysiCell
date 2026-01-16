#include "../core/PhysiCell.h"
#include "../modules/PhysiCell_standard_modules.h" 

using namespace BioFVM; 
using namespace PhysiCell;

void create_cell_types( void );
void setup_microenvironment( void ); 
void setup_tissue( void ); 
void tumor_cell_phenotype( Cell* pCell, Phenotype& phenotype, double dt );
std::vector<std::string> my_coloring_function( Cell* pCell );

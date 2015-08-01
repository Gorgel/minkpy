//Resample an array to a new resolution
void resampleArray(void *inArrayv, int inNumCells, int outNumCells, void *outArrayv)
{
        double *inArray = (double*)inArrayv;
        double *outArray = (double*)outArrayv;

        double inCellSize = 1.;
        float outCellSize = (double)inNumCells/(double)outNumCells*inCellSize;

	/* set default values */
	par.sigma=2; par.seed=200294;
	par.dim=intvector(1,3); par.dim[1]=64; par.dim[2]=32; par.dim[3]=16; 
	empty(&par.inname); empty(&par.outname);
	empty(&par.colname); empty(&par.grsname);
	par.length=1; par.nongauss=0; par.lo=-4; par.hi=4; par.bins=100; 
	par.med=0; par.pixel=1; par.normal=1; par.time=0;
	par.cstyle  =0; par.floattyp=0; par.complement = 0; 


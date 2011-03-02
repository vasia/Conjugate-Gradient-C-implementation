/*read random systems files*/

#include <stdio.h>
#include <stdlib.h>

//reads .sys to find N(number of buses), gen(number of generators) and br(number of branches) and SBASE
void readSys(int *N, int *gen, int *br, int *SBASE, FILE* f){	
	char *str;
	int first, I, Ig, Ib; char c;
	long Boffset, Goffset, Broffset;	//offset bus data, generator data and branch data
	*N = 0; *gen = 0; *br = 0;
	str = malloc(150*sizeof(char));

	//scanning bus data
	fscanf(f, "%d%c%d", &first, &c, &(*SBASE));
	fgets(str, 60, f);
	fgets(str, 60, f);
	fgets(str, 60, f);
	
	Boffset = ftell(f);	//start of bus data

	while(1) {			//BUS DATA
		fscanf(f, "%d", &I);
		if(I==0) {
			fgets(str, 100, f);
			Goffset = ftell(f);	//start of generator data
			break;
		}
		fgets(str, 100, f);
		(*N)++;	//counting number of buses
	}

	//scanning generator data
	fseek(f, Goffset, SEEK_SET);	//start of generator data
	while(1) {			
		fscanf(f, "%d", &Ig);
		if(Ig==0) {
			fgets(str, 150, f);
			Broffset = ftell(f);	//start of branch data
			break;
		}
		fgets(str, 150, f);
		(*gen)++;
	}

	//scanning branch data
	fseek(f, Broffset, SEEK_SET);	//start of branch data
	while(1) {			
		fscanf(f, "%d", &Ib);
		if(Ib==0) {
			fgets(str, 150, f);
			break;
		}
		fgets(str, 150, f);
		(*br)++;
	}
	

	free(str);

}

//reads .ses to find M(total number of measurements), voltageCount(number of voltage measurements), pfCount(number of real/reactive flow measurements) 
//and pCount(number of real/reactive injection measurements)
void readSes(int *M, int *voltageCount, int *pfCount, int *pCount, FILE* f){
	int Iv, Irf, Iri, Iqi, Iqf;
	char *str;
	*M = 0; *voltageCount = 0; *pfCount = 0; *pCount = 0;
	str = malloc(150*sizeof(char));
	long Voffset, RFoffset, QFoffset, RIoffset, QIoffset;	//offset voltage measurements, real flow, reactive flow, real injection, reactive injection

	fgets(str, 100, f);
	fgets(str, 100, f); 
	fgets(str, 100, f); 
	fgets(str, 100, f); 
	fgets(str, 100, f);

	Voffset = ftell(f);	//start of voltage measurements
	while(1) {			
		fscanf(f, "%d", &Iv);
		if(Iv==0) {
			fgets(str, 100, f);
			fgets(str, 100, f);
			RFoffset = ftell(f);	//real flow offset
			break;
		}
		fgets(str, 150, f);
		(*voltageCount)++;	//count total number of measurements
	}

	while(1) {			//real flow measurements
		fscanf(f, "%d", &Irf);
		if(Irf==0) {
			fgets(str, 100, f);
			fgets(str, 100, f);
			QFoffset = ftell(f);	//reactive flow offset
			break;
		}
		fgets(str, 150, f);
		(*pfCount)++;	//count number of flow measurements
	}
	
	while(1) {			//reactive flow measurements
		fscanf(f, "%d", &Iqf);
		if(Iqf==0) {
			fgets(str, 100, f);
			fgets(str, 100, f);
			RIoffset = ftell(f);	//real injection offset
			break;
		}
		fgets(str, 150, f);
	}

	while(1) {			//real injection measurements
		fscanf(f, "%d", &Iri);
		if(Iri==0) {
			fgets(str, 100, f);
			fgets(str, 100, f);
			QIoffset = ftell(f);	//reactive injection offset
			break;
		}
		fgets(str, 150, f);
		(*pCount)++;	//count number of injection measurements
	}
	
	while(1) {			//reactive injection measurements
		fscanf(f, "%d", &Iqi);
		if(Iqi==0) {
			fgets(str, 100, f);
			fgets(str, 100, f);
			break;
		}
		fgets(str, 150, f);
	}

	(*M)=*voltageCount+2*(*pfCount)+2*(*pCount);
	free(str);
} 



//reads .sys and stores bus, generator & branch data
void storeSys(struct sys *Sys, int N, int gen, int br, int *swing, FILE* f){
	double GL, BL, VM, VA, PG, QG, QT, QB, VS, MBASE, ZR, ZX, RT, XT, GTAP, RMPCT, PT, PB, RATEA, RATEB, RATEC;
	double GI, GJ, BI, BJ;
	char *str;
	int i, ZONE, IREG, fromtmp, totmp, gentmp, CKT, index;
	str = malloc(150*sizeof(char));
	fgets(str, 150, f);
	fgets(str, 150, f);
	fgets(str, 150, f);

	//scan bus data
	for(i=0; i<N; i++){
		fscanf(f, "%d %d %lf %lf %lf %lf", &(Sys->RBuses[i]), &(Sys->IDE[i]), &(Sys->PL[i]), &(Sys->QL[i]), &GL, &BL);
		
		//find swing virtual bus number
		if(Sys->IDE[i]==3) *swing=i;

		fscanf(f, "%d %lf %lf", &(Sys->IA[i]), &VM, &VA);
		name(&(Sys->Names[i+12]), f);
		fscanf(f, "%lf %d", &(Sys->BASKV[i]), &ZONE);
	}
	sort(Sys->RBuses, N, Sys->busIndex);//sort RBuses

	fgets(str, 100, f);	
	fgets(str, 100, f);	//move to generator data
	//scan generator data
	for(i=0; i<gen; i++){
		fscanf(f, "%d %d %lf%lf%lf%lf%lf%d%lf", &gentmp, &(Sys->ID[i]), &PG, &QG, &QT, &QB, &VS, &IREG, &MBASE );
		fscanf(f, "%lf%lf%lf%lf%lf%d%lf%lf%lf", &ZR, &ZX, &RT, &XT, &GTAP, &(Sys->STAT[i]), &RMPCT, &PT, &PB);
		//mapping
		index = BinSearch(Sys->RBuses,0, N, gentmp);
		Sys->Gbus[i] = index+1;
	}

	fgets(str, 150, f);	
	fgets(str, 150, f);	//move to branch data
	//scan branch data
	for(i=0; i<br; i++){
		fscanf(f, "%d %d %d %lf %lf %lf %lf %lf %lf", &fromtmp, &totmp, &CKT, &(Sys->R[i]), &(Sys->X[i]), &(Sys->BS[i]), &RATEA, &RATEB, &RATEC);
		fscanf(f, "%lf %lf %lf %lf %lf %lf %d", &(Sys->ratio[i]), &(Sys->angle[i]), &GI, &BI, &GJ, &BJ, &(Sys->st[i]));		

		//mapping
		index = BinSearch(Sys->RBuses,0, N, fromtmp);
		Sys->from[i] = index+1;
		index = BinSearch(Sys->RBuses,0, N, totmp);
		Sys->to[i] = index+1;

	}

	free(str);
}



//reads .ses and stores measurements data 
void storeSes(struct ses *Ses, struct Rvalues *rv, int *RBuses, int N, int voltageCount, int pfCount, int pCount, FILE* f){
	int Voltmp, index, i, fromtmp, totmp, CKT, ST, RTU, injtmp;
	char *str;
	str = malloc(150*sizeof(char));

	fgets(str, 100, f);
	fgets(str, 100, f); 
	fgets(str, 100, f); 
	fgets(str, 100, f); 
	fgets(str, 100, f);
	
	//store voltage measurements
	for(i=0; i<voltageCount; i++){
		fscanf(f, "%d %f %f %d %d %lf %lf", &Voltmp, &(rv->SNM[i]), &(rv->FS[i]), &ST, &RTU, &(rv->values[i]), &(rv->dev[i])); 
		fgets(str, 100, f);
	//mapping
		index = BinSearch(RBuses,0, N, Voltmp);
		Ses->busNumVol[i] = index+1;
	}
	
	fgets(str, 100, f); 
	fgets(str, 100, f);

	//store real flow measurements	
	int ii;
	for(i=0; i<pfCount; i++){
		ii = i+voltageCount;
		fscanf(f, "%d %d %d %f %f %d %d %lf %lf", &fromtmp, &totmp, &CKT, &(rv->SNM[ii]), &(rv->FS[ii]), &ST, &RTU, &(rv->values[ii]), &(rv->dev[ii])); 
		fgets(str, 100, f);
	//mapping
		index = BinSearch(RBuses,0, N, fromtmp);
		Ses->rfFrom[i] = index+1;
		index = BinSearch(RBuses,0, N, totmp);
		Ses->rfTo[i] = index+1;
	}
	
	fgets(str, 100, f); 
	fgets(str, 100, f);

	
	//store reactive flow measurements

	for(i=0; i<pfCount; i++){
		ii = i+voltageCount+pfCount;
		fscanf(f, "%d %d %d %f %f %d %d%lf %lf", &fromtmp, &totmp, &CKT, &(rv->SNM[ii]), &(rv->FS[ii]), &ST, &RTU, &(rv->values[ii]), &(rv->dev[ii])); 
	}

	fgets(str, 100, f);
	fgets(str, 100, f); 
	fgets(str, 100, f);

	
	//store real injection measurements
	for(i=0; i<pCount; i++){
		ii = i+voltageCount+2*pfCount;
		fscanf(f, "%d %f %f %d %d %lf %lf", &injtmp, &(rv->SNM[ii]), &(rv->FS[ii]), &ST, &RTU, &(rv->values[ii]), &(rv->dev[ii])); 
		fgets(str, 100, f);
	//mapping
		index = BinSearch(RBuses,0, N, injtmp);
		Ses->rBusNum[i] = index+1;
	}

	fgets(str, 100, f);
	fgets(str, 100, f); 

	//store reactive injection measurements
	for(i=0; i<pCount; i++){
		ii = i+voltageCount+2*pfCount+pCount;
		fscanf(f, "%d %f %f %d %d %lf %lf", &injtmp, &(rv->SNM[ii]), &(rv->FS[ii]), &ST, &RTU, &(rv->values[ii]), &(rv->dev[ii])); 
		fgets(str, 100, f);
	}

	free(str);
}


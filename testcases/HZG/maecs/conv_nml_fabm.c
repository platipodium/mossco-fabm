/*****************************************************/
/*  Simple namelist (nml) Parser for FABM and F2003  */
/*                                                   */
/*   Kai W. Wirtz  (HZG)                   10/06/13  */
/*****************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<string.h>
#define MI 9   		/* do not change these array sizes  ...   */
#define MAXP 59 	/*   unless you create highly complicated */
#define NAML 59 	/*   model and directory environments     */
/* ------------------------------------------------------------------------------ */
/*       configure-switches for the models; you MAY edit here:                    */ 
#define NewModF90 0	/* 1: creates an separate initialisation file; DO NOT USE */
#define TYPES 1		/* 1: creates an separate type file; obsolete in new FABM */
#define DIAG_RHS 0	/* 1: output of right hand sides as diagnostic variables  */
#define SED 0           /* 1: arranges extra calls for sediment models            */
#define SCALEFAC 0      /* 1: transorms all rate parameters already during get(.. */
                        /* 0: insert UNIT into _SET_ODE_(.. */
#define AUX 1		/* 1: adds derived parameters to the model structure      */
/* ------------------------------------------------------------------------------ */

char *fil(char *,char *, int ),*insul(char *);

main(int argn, char **argv )
{
/* ----------------------------------------------------------------------- */
/*  configuration variables for the model structure; you SHOULD edit here: */
/* ----------------------------------------------------------------------- */
char modname[NAML] ="maecs";/* model name; this string will be also used for building the nml-filenames */

char typet_nm[NAML]="hzg_";/* model-tree name for the type definition; FABM-convention includes the    */
				/* "institutional" origin, corresponding to folders within fabm/src/model/... */

char elements[7]   = "CNPS";	/* elements kept explicitly. other variants : "N", "CNP", or 
				/* "CNPSF" with "S" for Si/silicon or "F" for Fe/iron */
char longelemnam[7][NAML]={"carbon","nitrogen","phosphorus","silicate"};				

char nmlname0[MI][NAML] = {"switch","init","pars","graz","env","omex"}; /* list of all nml-files tags of the fabm-model,"ctl" */
				/* e.g., the tag "par" creates the full filename "modname_par.nml"  */
				/* MUST include an "init" for the definition/initialisation of state variables */
				/* MAY include a "switch" nml for all model switches to build modular hierarchies ; at the start of the list */
				
int nir=6;  			/* number of actual nml files; may be smaller than length of nmlname0 */
int bdim=1;  			/* number of spatial boxes in 0D */
char lstname0[MI][NAML] = {"deps","diags","aux"}; /* names of the files for 1) external forcing and 2) diagnostics and, optional (AUX=1), 3) derived parameters */
		/* order is fixed; ""deps" creates the full filename "modname_deps.lst"  */
				
char FabmDepVarName[NAML]= {"standard_variables%"}; /* prefix of the variable names in the FABM driver host */ 	
				/*  standard_variables\%*/
				
char dirn_nml[NAML]	= "";	/* directory where all the nml files  reside ./*/
char dirn_f90[3*NAML] 	= "/home/wirtz/mossco/fabm/src/models/hzg/maecs/";	/* directory where all the input sources (model.F90,...) reside ./
//char dirn_f90[3*NAML] = "/home/onur/opt/src/fabm-code/src/models/hzg/maecs/";	/* directory where all the input sources (model.F90,...) reside ./*/
char indent0[NAML] 	= "";	/* indentation in declaration part */ 
char init_incl[2*NAML] 	= "maecs_incl.lst"; /* file included in init routine*/
char init_varincl[3*NAML]= "call maecs_init_stoichvars(self)"; /* line included in init routine*/
//char vstructn[NAML]= "env";	/* name of the major variable structure */ 
char vstructn[NAML]	= "env";/* name of the major variable structure */ 
char vstructc[4]   	= "E";	/* elements major variable structure
          A: all state variables T: traits  E: environmental forcing N: nutrients */ 
char env_add[3*NAML] 	= "";	/* additional variable as member of the env structure (e.g., non-mass,non-trait MAECS variables*/ 
/* dependencies in switches */
char swi1[7][NAML]	= {"BGC2DDiagOn","Budget2DDiagOn","PhysiolDiagOn","BGC0DDiagOn","BGC0DDiagOn","-"};
char swi2[7][NAML]	= {"BGC0DDiagOn","Budget0DDiagOn","PhotoacclimOn","PhosphorusOn","BioOxyOn",""};

/* ------------------------------------------------------------------- */
/*            do not edit below ...                                    */
/* ------------------------------------------------------------------- */
int eoi,pi,pj,pjs,sws,ni,nj,d,np,nl,ss,i,pv,nmli=0,out=1,tti,nls,swi[MI][MAXP],
  setvel[MI][MAXP],setvel_n[MI][MAXP],trait[MAXP],found_rhs[MAXP],found_rate[MI][MAXP],pc;
unsigned long dind;
char line[256],line1[256],line2[256],*lr,c,*cp,*cp1,*cp2,*cp3,tabs[2][6]={"    ",""};
FILE *sp,*sp1,*sp3,*sp2,*sp4,*spv[3],*spt;
char keys[5][4]={"RHS","ODE","GET"};
// char traitpre[5]="phy%";/* structure name for functional group variables */
char insname[NAML]="";
char pref[NAML],pnam[NAML],ptyp[NAML],pvals[NAML],outn[256],outn2[256],ttn[NAML],fnam[NAML];
char nmlname[MI][NAML],nmlfname[MI][NAML];//"switch",
char nmlt[3][NAML] = {"par","init","switch"}; // 3 nml groups: parameters, init=states, logical switches
char numt[4][NAML] = {"real(rk)","logical","integer","character(len=64)"};
char yesno[2][NAML] = {"false","true"};
char tdb[2][4]={"","_sf"},tcb[2]={',','\n'};
char parname[MI][MAXP][NAML],snameshort[MAXP][NAML],snameshort2[MAXP][NAML],
  snameshort3[MAXP][NAML],parvals[MI][MAXP][NAML],pcom[MI][MAXP][3*NAML],texsymb[MI][MAXP][2*NAML],
  punit[MI][MAXP][NAML],sl[MAXP],setveln[MI][MAXP][NAML],swin[MI][MAXP][NAML],idname[MAXP][NAML],
  strvari[MAXP],partypen[MI][MAXP][NAML],pmapstring[MI][MAXP][2*NAML],*ptr,velem;
//char extvar[5][2][6]={{"I_0","par"},{"par","par"},{"",""}};
char sc[99], tmpstr[NAML];
int nmltype[MI],partype[MI][MAXP],nump[MI],ChemSpec[MAXP], ni0,nis,nelements;

// cp maecs_do_gen.F90 maecs/maecs_do.F90; cp maecs_types_gen.F90 maecs/maecs_types.F90; cp maecs_gen.F90 maecs/maecs.F90
 
/* printf("argn=%d\t",argn);printf("argv=%s\t",argv[1]);
if (argn>=2)   strcpy(simfile,argv[1]); else  exit(0);
if (argn>=3)   strcpy(cas,argv[2]);printf("AS=%s\n",cas);*/
nelements = strlen(elements);
tti=-1;
// ----------------------------------------------------- 
//  read all info from comment header in nml files
// ----------------------------------------------------- 
for(ni=0,nis=-1;ni<nir;ni++)
  {
  if(out) printf("\nstart %d\n",ni);
// ----------------------------------------------------- 
//  create names of nml files
// ----------------------------------------------------- 
  sprintf(nmlfname[ni],"%s%s_%s.nml",dirn_nml,modname,nmlname0[ni]);
  sprintf(nmlname[ni],"%s_%s",modname,nmlname0[ni]);
// ----------------------------------------------------- 
//  detect type of nml files
// ----------------------------------------------------- 
  nmltype[ni]=0;  /* standard par list */
  if(strstr(nmlname0[ni],"init")!=NULL) nmltype[ni]=1,ni0=ni; /* list of state variables */
  if(strcmp(nmlname0[ni],"switch")==0) nmltype[ni]=2,nis=ni; /* list of logical switches */
    
  if(out) printf("extracting namelist from %s group-type %s \n",nmlfname[ni],nmlt[nmltype[ni]]);

// ----------------------------------------------------- 
//  open nml files for reading
// ----------------------------------------------------- 
  sp=fopen(nmlfname[ni],"r"); 
  if(sp==NULL) {printf("Problem while opening %s !\n",nmlfname[ni]),exit(0);}
  
// ----  reads 3 initial comment lines 
  for(i=0;i<4;i++) lr=fgets(line,256,sp); 

  lr=(char *)1; pi=0; eoi=0;
  while(lr!=NULL && pi<MAXP && eoi==0)
    {
// ----  retrieve info line per line 
    sscanf(line,"%s %s %s",pref,pnam,ptyp);
// ----  check for units
    cp=strchr(line,'='); 
    if(cp!=NULL)
      {
      cp1=cp;  cp=cp+2; c=*cp; 
      while( (c==' ' || c=='\t') && cp-cp1<254)
        c=*(++cp);
      if((cp1=strstr(cp,"]"))>=0)
         strncpy(punit[ni][pi],cp,cp1-cp);
//       else if((cp1=strstr(cp,"\n"))>=0)
//          strncpy(punit[ni][pi],cp,cp1-cp);
      else  
         strcpy(punit[ni][pi],cp);   
      }
    else
      strcpy(punit[ni][pi],"");
    
//  ----  check for switch-dependency
    swi[ni][pi]=-1;
    cp=strchr(line,'#'); 
    if(cp!=NULL && nis>=0)
      {
      cp1=strchr(line,'\n'); 
      strncpy(swin[ni][pi],cp+1,cp1-cp-1);
    
      if (out) printf("*** switch-dependency %s %s...\n",pnam,swin[ni][pi]);
      for(pj=0;pj<nump[nis];pj++)
        if(strcmp(swin[ni][pi],parname[nis][pj])==0) swi[ni][pi]=pj, pj=nump[nis]+3;  
      if(pj<nump[nis]+3)
        {
	printf("\n%s not found in switch-namelist !!\n\n",swin[ni][pi]);
        for(pj=0;pj<nump[nis];pj++)
	  printf("%s %d\t",parname[nis][pj],strcmp(swin[ni][pi],parname[nis][pj]));
        printf("\n\n");
	}
      }
      
//  ----  check for trait transporter (e.g., phyC)
    if(strstr(line,"!TT") && ni==ni0)  
       tti=pi;    
//  ----  check for env variables (not in maecs structure)
    if(strstr(line,"EE") && ni==ni0) 
      {
      strcat(env_add,pnam),
      cp=rindex(env_add,'_');*cp='\0';
      strcat(env_add,", ");
      }
      
//  ----  check for tech symbols needed for doxygen
    strcpy(texsymb[ni][pi],""); 
    cp1=strchr(line,'$'); 
    cp =strrchr(line,'$'); 
    
    if(cp1!=NULL)
      {
       if(cp>cp1+1)
         strncpy(texsymb[ni][pi],cp1+1,cp-cp1-1);
     }
    else
      {
      //if no symbol was provided, texsymb[ni][pi] should take: '\mathrm{pnam}
      memset(&tmpstr[0], 0, sizeof(tmpstr));
      strcpy(tmpstr,"\\mathrm{");
      strcat(tmpstr,pnam);
      strcat(tmpstr,"}");
      strcpy(texsymb[ni][pi],tmpstr);
      }     
    
    lr=fgets(line,256,sp);
    cp=strchr(line,'!'); cp1=cp;
    cp=cp+2; c=*cp; 
// ----- search start of par name string 
    while( (c==' ' || c=='\t') && cp-cp1<254)
       c=*(++cp);
    if((cp1=strstr(cp,"initial "))>=cp)
        cp+=strlen("initial "); 
 
    if((cp1=strstr(cp,"\n"))>=0)
       strncpy(pcom[ni][pi],cp,cp1-cp);
    else
       strcpy(pcom[ni][pi],cp);
    
    //  ----  check for settling velocity
    setvel[ni][pi]=-1;
//    printf("%s %ld\n",cp,strstr(line,"velocity"));
    if((cp=strstr(pcom[ni][pi],"velocity")) !=NULL )
      {
      cp+=strlen("velocity")+1;
//      printf("compare velocity identities:\t %s ...\n",cp);
      for(pj=0;pj<nump[ni0];pj++)
        {
	strcpy(line2,parname[ni0][pj]);strcat(line2," ");
        strcat(line2,pcom[ni0][pj]);
//         printf("\t ...with %s ...\n",line2);
	ptr = strtok(line2," ");
        while(ptr != NULL)
	  {
//	  printf("check state word %s ...\n",ptr);
          if(strstr(cp,ptr) !=NULL)
	    {
	  printf("***match %s %s ... linked to state %d\n",cp,ptr,pj);
	    setvel[ni][pi] = pj;
	    setvel[ni0][pj]= pi; setvel_n[ni0][pj]=ni;
//	    pj=nump[ni0]+3;  
   	    }
	  ptr = strtok(NULL," ");
	  }
	}  
      if(pj<nump[nis]+3) printf("\n %s not found in state-namelst !!\n\n",cp);
      if (out) printf("*** settling velocity %s linked to state %d\t %s...\n",pnam,setvel[ni][pi],parname[ni0][setvel[ni][pi]]);
      }   
    
    strcpy(parname[ni][pi],pnam);    
    
    if(strstr(ptyp,"bool")) partype[ni][pi]=1;
    else partype[ni][pi]=0;
    if(strstr(ptyp,"integer")) partype[ni][pi]=2;
    
 
    if(out) printf("%d %s %d unit=%s\t # %s\n",pi,parname[ni][pi],partype[ni][pi],punit[ni][pi],pcom[ni][pi]);
    pi++;
    lr=fgets(line,256,sp);
    if(strstr(line,"!--")) eoi=1; 
    }
/*   
for(d=0;d<nmi;d++)   
  fprintf(sp,"\n\n open(namlst,file='%s',status='old')\n read(namlst,nml=%s',err=%d)\n",nlname,fn);  */ 
  while(lr!=NULL && line[0]!='&')
     lr=fgets(line,256,sp); 
  sscanf(&line[1],"%s",pnam);
  if(strcmp(pnam,nmlname[ni])!=0)
    {
    printf("\n ***\n The input namelist name \"%s\" does not match \"%s\" !!\n ***\n",pnam,nmlname[ni]);
    printf("\n  please edit, and re-run the parser ...!!\n\n");
    exit(0);
    }
  else if(out) printf("namelist %s = %s\n",pnam,nmlname[ni]);

// ----------------------------------------------------- 
//  read parameter values as strings, not used yet
// ----------------------------------------------------- 
  lr=fgets(line,256,sp);
  pv=0;  eoi=0;
  while(lr!=NULL && eoi==0)
    {
    sscanf(line,"%s %s %s",pnam,pref,pvals);
//   printf("found %s %s\n",pnam,pvals);
    for(pj=0;pj<pi;pj++)
      if(strcmp(pnam,parname[ni][pj])==0)
         {
         if((cp1=strstr(pvals,","))!=NULL)
            strncpy(parvals[ni][pj],pvals,cp1-pvals);
	 else
            strcpy(parvals[ni][pj],pvals);
	 if(out) printf("found %s for %d\tval=%s\n",pnam,pj,parvals[ni][pj]);
	 pj=pi;
         pv++;
	 }
    lr=fgets(line,256,sp);
    if(strstr(line,"/")!=NULL) eoi=1; 
    }
// ----------------------------------------------------- 
//  consistency check (header, value list)
// ----------------------------------------------------- 
  if(pv!=pi)
     printf("number of values %d doesn't match parameter list %d\n",pv,pi);
  else if(out)
     printf("found %d parameter values\n",pv);
    
  nump[ni]=pi;  
  fclose(sp);
  }
// ----------------------------------------------------- 
//  detect trait/fractional variables
// ----------------------------------------------------- 
ni=ni0;  
for(pj=0;pj<nump[ni];pj++)
    {
    strcpy(pnam,parname[ni][pj]);
    d=0,trait[pj]=0; 
    if(strstr(pnam,"frac_")) d=5,trait[pj]=1;  
    if(strstr(pnam,"size")) trait[pj]=1;
    cp=rindex(pnam,'_');
    *cp='\0';
// all state variables not ending with an element symbol are assumed to be traits 
    if(strchr(elements,pnam[strlen(pnam)])==NULL) 
      trait[pj]=1; 
    if(trait[pj] && out) printf("%s identified as trait\n",pnam);
    strcpy(snameshort[pj],pnam);
// ----------------------------------------------------- 
//  remove prefix in trait/fractional variable name 
// ----------------------------------------------------- 
    strcpy(snameshort2[pj],&pnam[d]);
    strcpy(snameshort3[pj],&pnam[d]);
    found_rhs[pj]=0;
    }
    
// ----------------------------------------------------- 
//   check for trait transporter (e.g., phyC)
// ----------------------------------------------------- 
if(tti>=0)
  strcpy(ttn,parname[ni][tti]);
else
  strcpy(ttn,"Biomass");


// ----------------------------------------------------- 
//  read all info from comment header in nml files
// ----------------------------------------------------- 
for(nls=0,ni=nir;nls<2+AUX;nls++,ni++)
  {
// ----------------------------------------------------- 
//  create filename 
// ----------------------------------------------------- 
  sprintf(fnam,"%s%s_%s.lst",dirn_nml,modname,lstname0[nls]);
  if(out) printf("extracting %s list from %s ...\n",lstname0[nls],fnam);
// ----------------------------------------------------- 
//  open lst files for reading
// ----------------------------------------------------- 
  sp=fopen(fnam,"r"); 
  if(sp==NULL) {printf("Problem while opening %s !\n",fnam),exit(0);}
  
// ----  reads 3 initial comment lines 
  for(i=0;i<4;i++) lr=fgets(line,256,sp); 
  lr=(char *)1; pi=0; eoi=0;
  while(lr!=NULL && pi<MAXP && eoi==0)
    {
// ----  retrieve info line per line 
//    sscanf(line,"%s %s %s %s %[ a-z:-#A-Z]\n",parname[ni][pi],punit[ni][pi],
    sscanf(line,"%s %s %s %s %s\n",parname[ni][pi],punit[ni][pi],
	   partypen[ni][pi],pmapstring[ni][pi],pcom[ni][pi]);
    if(out) printf("%d %s unit=%s\t\tpc=%s\n",pi,parname[ni][pi],punit[ni][pi],pcom[ni][pi]);
//  ----  check for switch-dependency
    swi[ni][pi]=-1;
    cp=strchr(pcom[ni][pi],'#'); 
    if(cp!=NULL && nis>=0)
      {
      strcpy(swin[ni][pi],cp+1);
      if(out) printf("*** switch-dependency %s:  %s (%ls, %s)\n",parname[ni][pi],swin[ni][pi],cp,cp);
      for(pj=0;pj<nump[nis];pj++)
        if(strcmp(swin[ni][pi],parname[nis][pj])==0) swi[ni][pi]=pj, pj=nump[nis]+3;  
      if(pj<nump[nis]+3) printf("\nnot found in switch-namelist !!\n\n");
      *cp='\0';
      }
    lr=fgets(line,256,sp);
    pi++;  
    }
//    fprintf(sp,"adding %d rhs to diagnostics\n",ni0);
  if(DIAG_RHS && strcmp(lstname0[nls],"diag")>=0)
    for(pj=0;pj<nump[ni0];pj++,pi++)
      {
      sprintf(parname[ni][pi],"rhs_%s",snameshort2[pj]);
      sprintf(punit[ni][pi],"%s/d",punit[ni0][pj]);
      strcpy(partypen[ni][pi],"step_integrated");
//      strcpy(partypen[ni][pi],"last");
      sprintf(pmapstring[ni][pi],"rhsv%c%s",'%',snameshort2[pj]);
      sprintf(pcom[ni][pi],"RHS of %s",pcom[ni0][pj]);
      }
  nump[ni]=pi;   
  }
  //  ----  check for rates
for(ni=0;ni<nir+2;ni++)     
   for(pj=0;pj<nump[ni];pj++)
     {
     found_rate[ni][pj]=-1;
     if(ni!=ni0 && ni!=nis)
       if(strstr(pcom[ni][pj],"rate") !=NULL || strstr(pcom[ni][pj],"RHS") !=NULL|| strstr(pcom[ni][pj],"velocity") !=NULL || parname[ni][pj][0]=='r' )
        {
        if (out) printf("*** found rate %s ...\n",parname[ni][pj]);
        found_rate[ni][pj]=1;
        }
     }
 printf("\n ------------------------\n %s \t%d\n\n",elements,nelements);

// ---------------------------------------------------------------- 
// ---------------------------------------------------------------- 
//                         output
// ---------------------------------------------------------------- 
// ---------------------------------------------------------------- 
sprintf(outn,"%s_gen.F90\0",modname); //initvar_
sp=fopen(outn,"w"); 
if(sp==NULL) {printf("Problem while opening %s !\n",outn),exit(0);}
if(out || 1) printf("writing definitions in %s...\n",outn);

// read existing code before listed initialisation  
sprintf(insname,"%s%s.F90\0",dirn_f90,modname);
printf("reading from %s ...\n",insname);
sp1=fopen(insname,"r"); 
if(sp1==NULL) {printf("Problem while opening %s !\n",insname),exit(0);}
// if(strlen(insname)>1)  {
//sprintf(outn,"~/maecs/FABM/fabm-kse/src/models/hzg/%s/%s",modname,insname);
//fprintf(sp,"! code included from %s\n",insname);
lr=(char *)1;
cp1=strstr(modname,"#SP#");
while(lr!=NULL && cp1==NULL)
  {
  lr=fgets(line,256,sp1); 
  fputs(line,sp);
  cp1=strstr(line,"#SP0");
  }
  
if(TYPES)
  {
  sprintf(insname,"%s%s_types.F90\0",dirn_f90,modname);
  printf("reading from %s ...\n",insname);
  sp2=fopen(insname,"r"); 
  if(sp2==NULL) {printf("Problem while opening %s !\n",insname),exit(0);}
  sprintf(outn,"%s_types_gen.F90\0",modname);
  sp3=fopen(outn,"w"); 
  if(sp==NULL) {printf("Problem while opening %s !\n",outn),exit(0);}
  if(out || 1) printf("writing type definitions in %s...\n",outn);
 
  lr=(char *)1;
  cp1=strstr(modname,"#SP#");
  while(lr!=NULL && cp1==NULL)
    {
    lr=fgets(line,256,sp2); 
    fputs(line,sp3);
    cp1=strstr(line,"#SP0");
    }
  fprintf(sp3,"!!----------------------------------------------------------------------\n");
  fprintf(sp3,"!! this code is generated by a parser  (conv_nml_fabm.c by kai wirtz)\n");
  fprintf(sp3,"!!----------------------------------------------------------------------\n");
  spt=sp3;  
  }
  else
    spt=sp;
fprintf(sp,"!!----------------------------------------------------------------------\n");
fprintf(sp,"!! this code is generated by a parser  (conv_nml_fabm.c by kai wirtz)\n");
fprintf(sp,"!!----------------------------------------------------------------------\n");
/* fprintf(sp,"#include \"fabm_driver.h\"\n#ifdef _FABM_F2003_\n");
fprintf(sp,"module fabm_%s\n!\n",modname);
fprintf(sp,"! !USES:\n  use fabm_types\n  use fabm_driver\n");
if(TYPES) fprintf(sp,"  use %s_types\n  private\n",modname);
fprintf(sp,"  implicit none\n\n !public flags for communication with fabm-driver/time-loop\n");
fprintf(sp,"  character(len=80),   public      ::  trait_timestr\n");
fprintf(sp,"  integer,             public      ::  ShowTraitRGR\n  private\n");
*/

fprintf(sp,"! --- HZG model types\n");

// ---------------------------------------------------------------- 
//    define the major variable structure of the model 
fprintf(spt,"! standard fabm model types\n",indent0);

if(TYPES)
  sprintf(insname,"%s_base_model\n",modname);
else
  sprintf(insname,"%s%s\n",typet_nm,modname);
//if(strcmp(modname,"maecs")!=0)
//  fprintf(spt,"%stype,extends(type_base_model),public :: type_%s",indent0,insname);
if (bdim<=1)
  fprintf(spt,"%stype (type_state_variable_id)        :: ",indent0);
else
  fprintf(spt,"%stype (type_state_variable_id), dimension(%d) :: ",indent0,bdim);
 
for(ni=ni0,pj=0;pj<nump[ni];pj++)
  fprintf(spt,"id_%s%c",snameshort2[pj],tcb[pj==nump[ni]-1]);

for(ni=nir,pj=0;pj<nump[ni];pj++)
  fprintf(spt,"%stype (type_%s_id)            :: id_%s\n",indent0,partypen[ni][pj],parname[ni][pj]);

fprintf(spt,"%stype (type_diagnostic_variable_id)   ::",indent0);
for(ni=nir+1,pj=0;pj<nump[ni];pj++)
  fprintf(spt," id_%s%c",parname[ni][pj],tcb[pj==nump[ni]-1]);


for(ni=0;ni<nir;ni++)
 if(ni!=nis )//&& ni!=ni0
   {
   fprintf(spt,"%sreal(rk) :: ",indent0);
   for(pj=0;pj<nump[ni];pj++)
     fprintf(spt," %s%c",parname[ni][pj],tcb[pj==nump[ni]-1]);
   }
   
// derived parameters 
if(AUX)
  {
  fprintf(spt,"%sreal(rk) :: ",indent0);
  for(ni=nir+2,pj=0;pj<nump[ni];pj++)
    fprintf(spt," %s%c",parname[ni][pj],tcb[pj==nump[ni]-1]);
  }

if(nis>=0)
  {
  fprintf(spt,"%slogical  :: ",indent0);
  for(ni=nis,pj=0,pi=0;pj<nump[ni];pj++)
   if(partype[ni][pj]==1)
     fprintf(spt," %s%c",parname[ni][pj],tcb[pj==nump[ni]-1]);
   else
     pi++;
  if(pi>0)
    {
    fprintf(spt,"%sinteger  :: ",indent0);
    for(ni=nis,pj=0,pjs=0;pj<nump[ni];pj++)
    if(partype[ni][pj]==2)
     fprintf(spt," %s%c",parname[ni][pj],tcb[(pjs++)==pi-1]);
    }
  }
    
if(TYPES==0)
  {
  fprintf(sp,"\n%scontains\n!   Model procedures\n",indent0);
  fprintf(sp,"%sprocedure :: initialize\n%sprocedure :: do\n",indent0,indent0);
  if(bdim>1) fprintf(sp,"%sprocedure :: mixing\n%sprocedure :: errfunc\n",indent0,indent0);


//fprintf(sp,"  procedure :: do_ppdd\n%sprocedure :: get_light_extinction\n%sprocedure :: get_conserved_quantities\n%sprocedure :: out_RGR_trait\n\n");
  }
// else
//  sp=sp3;
fprintf(spt,"%send type type_%s\n!\n",indent0,insname);

// ---------------------------------------------------------------- 
//    define the major variable structure of the model 
if(strlen(vstructn) > 0)
  {
  fprintf(spt,"%stype type_%s_%s\n",indent0,modname,vstructn);
  if (strchr(vstructc,'A')!=NULL)
    {
//    fprintf(spt,"%s real(rk) :: ",indent0);
    fprintf(spt,"%s real(rk) :: %s ",indent0,env_add); 
    for(ni=ni0,pj=0;pj<nump[ni];pj++)
      fprintf(spt,"%s%c",snameshort2[pj],tcb[pj==nump[ni]-1]);
    }
  if (strchr(vstructc,'E')!=NULL) 
    {
    fprintf(spt,"%s real(rk) :: %s ",indent0,env_add); 
    for(ni=nir,pj=0;pj<nump[ni];pj++)
      if(strstr(parname[ni][pj],"tot")==NULL)
        fprintf(spt,"%s%c",parname[ni][pj],tcb[pj==nump[ni]-1]);
    }
 if (strchr(vstructc,'N')!=NULL)
    {
    fprintf(spt,"%s real(rk) :: ",indent0);
    for(ni=ni0,pj=0,pi=0;pj<nump[ni];pj++)
     if(strstr(snameshort2[pj],"nut")!=NULL)
       {
       if(pi>0) fprintf(spt,", ");  
       pi++; 
       fprintf(spt,"%s",snameshort2[pj]);  
       }
    }
  fprintf(spt,"%send type\n",indent0);
  }
  
  
// ---------------------------------------------------------------- 
//    define the RHS structure containibg all state variables
fprintf(spt,"%stype type_%s_rhs\n%s real(rk) :: ",indent0,modname,indent0);
for(ni=ni0,pj=0;pj<nump[ni];pj++)
  fprintf(spt,"%s%c",snameshort2[pj],tcb[pj==nump[ni]-1]);
fprintf(spt,"%send type\n! standard fabm model types\n",indent0);
  
//if(nelements>0) fprintf(spt,"%stype (type_conserved_quantity_id)    ::",indent0); 
//for(pj=0;pj<nelements;pj++)
//  fprintf(spt," id_tot%c%c",elements[pj],tcb[pj==nelements-1]);
    
if(TYPES==0)
  fprintf(sp,"! !PRIVATE DATA MEMBERS:\n real(rk), parameter :: secs_pr_day = 86400.0_rk\n!EOP\n");
fprintf(sp,"!!--------------------------------------------------------------------\n");
fprintf(sp,"\n%scontains\n",indent0,modname);
// this is old documentation style: fprintf(sp,"\n%scontains\n! !IROUTINE: Initialise the %s model\n",indent0,modname);

// documentation for doxygen
fprintf(sp,"!> @brief initializes the model\n",indent0);
fprintf(sp,"!! @details here the maecs namelists are read and assigned respectively in the model type (self),\n",indent0); 
fprintf(sp,"!! state & diagnostic variables are registered in FABM and dependencies are imported from FABM\n",indent0); 
fprintf(sp,"!>\n!> **Model parameters, descriptions and corresponding symbols used in formulas:**\n",indent0);
// /describepar list
fprintf(sp,"! initial values\n",indent0);
for(pj=0;pj<nump[ni];pj++) 
  {
  sprintf(line,"%s",texsymb[ni][pj]);
  fprintf(sp,"!> \\describepar{%s%s, %s%s, %s, %s %s}\n",insul(parname[ni][pj]),fil(sc,parname[ni][pj],13),line,fil(sc,line,22),pcom[ni][pj],parvals[ni][pj],punit[ni][pj]); 
  if(trait[pj])
    fprintf(sp,"!> \\describepar{%s%s, %s%s, trait x biomass}\n",snameshort2[pj],fil(sc,snameshort2[pj],13),line,texsymb[ni0][tti]);
  } 
fprintf(sp,"! other parameters\n",indent0);
for(ni=0;ni<nir;ni++)
 if(ni!=nis && ni!=ni0)
  for(pj=0;pj<nump[ni];pj++)
    {
    sprintf(line,"%s",texsymb[ni][pj]);
    sprintf(line2,"%s",pcom[ni][pj]);  
    fprintf(sp,"!> \\describepar{%s%s, %s%s, %s, %s %s}\n",insul(parname[ni][pj]),fil(sc,parname[ni][pj],13),line,fil(sc,line,22),line2,parvals[ni][pj],punit[ni][pj]);  
    }
    
//this contains old documentation style: fprintf(sp,"!\n! !INTERFACE:\n%ssubroutine initialize(self, configunit)\n",indent0);
fprintf(sp,"%ssubroutine initialize(self, configunit)\n\n",indent0);
// this is old documentation style: fprintf(sp,"! !DESCRIPTION:\n!  Here, the namelists are read and the variables exported\n");
// this is old documentation style: fprintf(sp,"!  by the model are registered with FABM.\n!\n! !INPUT PARAMETERS:\n");
fprintf(sp,"%sclass (type_%s%s), intent(inout), target :: self\n",indent0,typet_nm,modname);
fprintf(sp,"%sinteger,                  intent(in)            :: configunit\n",indent0);
fprintf(sp,"!\n! !LOCAL VARIABLES:\n");
fprintf(sp,"%sinteger    :: namlst=19\n",indent0);
/* write Initial values */

fprintf(sp,"!!------- Initial values of model %s ------- \n",modname);
ni=ni0;
if(out) printf("writing %d state variables from list %d ...\n",nump[ni],ni);


for(pj=0;pj<nump[ni];pj++)
  {
  fprintf(sp,"%s%s  :: %s%s! %s\n",indent0,numt[partype[ni][pj]],parname[ni][pj],fil(sc,parname[ni][pj],13),pcom[ni][pj]); /* ,parvals[ni][pj]*/
  if(trait[pj])
    fprintf(sp,"%s%s  :: %s  ! trait times biomass\n",indent0,numt[partype[ni][pj]],snameshort2[pj]);
//  if(out) printf("  %s\t:: %s=0.d0  ! %s\n",numt[partype[ni][pj]],parname[ni][pj],pcom[ni][pj]); /* ,parvals[ni][pj]*/
  } 
 
/* write all par values */
for(ni=0;ni<nir;ni++)
 if(ni!=nis && ni!=ni0)
  {
  fprintf(sp,"!!------- Parameters from nml-list %s ------- \n",nmlname[ni]);
  if(out) printf("writing %d parameters from list %d %s\n",nump[ni],ni,nmlname[ni]); 
  for(pj=0;pj<nump[ni];pj++)
    {
    sprintf(line,"%s",parname[ni][pj]);/*,parvals[ni][pj]*/
//    if(strlen(line)<5) strcat(line,"\t");
    fprintf(sp,"%s%s  :: %s%s! %s\n",indent0,numt[partype[ni][pj]],line,fil(sc,line,13),pcom[ni][pj]);    
//    if(out) printf("  %s\t:: %s  ! %s\n",numt[partype[ni][pj]],line,pcom[ni][pj]); 
    }
  }
 /* write switches  */
if((ni=nis)>=0)
  {
  fprintf(sp,"!!------- Switches for configuring model structure -------\n");
  if(out) printf("writing %d switches from list %d\n",nump[ni],ni);
  for(pj=0;pj<nump[ni];pj++)
    {
    sprintf(line,"%s",parname[ni][pj]);/*,parvals[ni][pj]*/
    fprintf(sp,"%s%s   :: %s%s! %s\n",indent0,numt[partype[ni][pj]],line,fil(sc,line,13),pcom[ni][pj]); 
//  if(out) printf("switch  %s\t:: %s  ! %s\n",numt[partype[ni][pj]],line,pcom[ni][pj]); 
    }
  }
if(bdim>1)  
   fprintf(sp,"%s   :: sbname,cbname\ninteger   :: ib\n",numt[2]);

for(ni=0;ni<nir;ni++)
  {
  fprintf(sp,"\n%snamelist /%s/ &\n%s  ",indent0,nmlname[ni],indent0);
  for(pj=0,pi=0;pj<nump[ni]-1;pj++)
    {
    fprintf(sp,"%s, ",parname[ni][pj]);
    pi+=strlen(parname[ni][pj]);
    if(pi>56 && pj<nump[ni]-1) fprintf(sp,"&\n%s  ",indent0), pi=0;
    }
  fprintf(sp,"%s\n",parname[ni][pj]);  
  }
fprintf(sp,"\n");  
/* Initialize all input variables  */
for(ni=0;ni<nir;ni++)
 if(ni!=nis)
  {
  for(pj=0;pj<nump[ni];pj++)
    {
    sprintf(line,"%s",parname[ni][pj]);/*,parvals[ni][pj]*/
    sprintf(line2,"%s",parvals[ni][pj]);/*,parvals[ni][pj]*/

    if(ni!=nis)
      {
      fprintf(sp,"%s%s%s= %s_rk",indent0,line,fil(sc,line,13),line2);
      fprintf(sp,"%s! %s\n",fil(sc,line2,16),punit[ni][pj]);
      }//fil(sc,line2,20) 
    else      
      fprintf(sp,"%s%s%s = %s\n",indent0,line,fil(sc,line,13),parvals[ni][pj]);  
    }
  }
  
fprintf(sp,"\n\n!--------- read namelists --------- \n");  
fprintf(sp,"write(0,*) ' read namelists ....'\n");  

for(ni=0;ni<nir;ni++) {
//  fprintf(sp,"\nopen(namlst,file='%s',status='old')",nmlfname[ni]);
//    fprintf(sp,"  read(configunit,nml=%s,err=%d,end=%d)\n",nmlname[ni],90+ni,90+ni+MI);
  fprintf(sp,"%sopen(namlst,file='%s',status='old')\n",indent0,nmlfname[ni]);
  fprintf(sp,"%sread(namlst,nml=%s,err=%d,end=%d)\n",indent0,nmlname[ni],90+ni,90+ni+MI);
}
fprintf(sp,"! Store parameter values in our own derived type\n");  
fprintf(sp,"! NB: all rates must be provided in values per day,\n");  
fprintf(sp,"! and are converted here to values per second.\n");  

if(nis>=0)
  {
  fprintf(sp,"\n!!------- logical parameters: switches  -------\n");
  ni=nis;
  
  for(pj=0;pj<nump[ni];pj++)
    {
    strcpy(line2, fil(sc,parname[ni][pj],13) );
    fprintf(sp,"%scall self%cget_parameter(self%c%s,%s '%s',%s default=%s)\n",indent0,'%','%',parname[ni][pj],line2,parname[ni][pj],line2,parname[ni][pj]); 
    }
  }
for(ni=0;ni<nir;ni++)
 if( ni!=nis)//ni!=ni0 &&
  {
  fprintf(sp,"\n!!------- model parameters from nml-list %s ------- \n",nmlname[ni]);
  for(pjs=-1;pjs<nump[nis];pjs++)
   for(pj=0,sws=0;pj<nump[ni];pj++)
    {
//    printf("%d %d\n",swi[ni][pj],pjs);
    if(swi[ni][pj]==pjs)
      {
      if(sws==0 && pjs>=0 && ni!=ni0)
         fprintf(sp,"if (self%c%s) then\n",'%',swin[ni][pj]), sws=1; 
      sprintf(line,"%s",parname[ni][pj]);/*,parvals[ni][pj]*/
      strcpy(line2, fil(sc,line,13) );
      strcat(line,line2);
      
//      if(out)printf("%sself%c%s= %s\n",tabs[(pjs<0)],'%',line,parname[ni][pj]);
      fprintf(sp,"%s%scall self%cget_parameter(self%c%s,'%s',%s default=%s", &
	   indent0,tabs[(pjs<0)],'%','%',line,parname[ni][pj],line2,parname[ni][pj]); 
      if(found_rate[ni][pj]==1 && SCALEFAC) fprintf(sp,", scale_factor=1.0_rk/secs_pr_day");
      fprintf(sp,")\n");
      }
    if( pj==nump[ni]-1 && sws==1) fprintf(sp,"end if\n");
    }
  }
  
 /* write derived parameters   */
if(AUX)
  {
  fprintf(sp,"\n!!------- derived parameters  ------- \n");
  for(ni=nir+2,pj=0;pj<nump[ni];pj++)
    fprintf(sp,"%sself%c%s%s = %s\n",indent0,'%',parname[ni][pj],fil(sc,parname[ni][pj],12),pmapstring[ni][pj]);
  }
  
ni=ni0;
printf("!------- Register state variables  ------- \n");
fprintf(sp,"\n!!------- Register state variables  ------- \n");

if (bdim>1)
  fprintf(sp,"do ib = 1, %d\n",bdim);

for(pjs=-1;pjs<nump[nis];pjs++)
  for(pj=0,sws=0;pj<nump[ni];pj++)
     {
     if(swi[ni][pj]==pjs)
      {
      sl[pj]= ' '; if(strlen(snameshort2[pj])<7) sl[pj]='\t';
        
      if(sws==0 && pjs>=0)
        fprintf(sp,"\nif (self%c%s) then\n",'%',swin[ni][pj]), sws=1;
      
      if(pjs>=0 &&out) printf("%s\t%s %s\t %d\n",swin[ni][pj],snameshort[pj],snameshort2[pj],trait[pj]);
     
      if(trait[pj])
        fprintf(sp,"%s%s = %s * %s  ! trait times biomass\n",tabs[(pjs<0)],snameshort2[pj],parname[ni][pj],ttn); /* ,parvals[ni][pj]*/

     if (bdim>1)
      {
//	fname = TRIM(out_dir) //'/'// TRIM(out_fn) // '.' // ext
      fprintf(sp,"write (sbname, \"(A,I1)\") '%s',ib \n",snameshort2[pj]);
      fprintf(sp,"cbname = '%s ' // sbname\n",pcom[ni][pj]);  // ib
      fprintf(sp,"%s%scall self%cregister_state_variable(self%cid_%s(ib),sbname,",indent0,tabs[(pjs<0)],'%','%',snameshort2[pj]);
      fprintf(sp,"'%s',cbname, &\n",punit[ni][pj]); 
      sprintf(idname[pj],"%s(ib)",snameshort2[pj]);
      }
     else
      {
      fprintf(sp,"%s%scall self%cregister_state_variable(self%cid_%s,%s'%s',",indent0,tabs[(pjs<0)],'%','%',snameshort2[pj],fil(sc,snameshort2[pj],6),snameshort2[pj]);
      fprintf(sp,"'%s','%s %s', &\n",punit[ni][pj],pcom[ni][pj],snameshort2[pj]); 
      sprintf(idname[pj],"%s",snameshort2[pj]);
      }
//      fprintf(sp,"%s\t\'%s\', \'%s\', &\n",tabs[(pjs<0)],punit[ni][pj],pcom[ni][pj]); 
      if(trait[pj])
        fprintf(sp,"%s   %s",tabs[(pjs<0)],snameshort2[pj]); 
      else
        fprintf(sp,"%s   %s",tabs[(pjs<0)],parname[ni][pj]); 
      
      ss=0;
      if(SED==0 && snameshort2[pj][0]=='o') ss=1; 
      if(strcmp(modname,"maecs")==0)
        {
        if(trait[pj]==1 || strncmp(snameshort2[pj],"phy",3)==0 || strncmp(snameshort2[pj],"zoo",3)==0)
          fprintf(sp,", minimum=_ZERO_, no_river_dilution=plankton_no_river_dilution "); 
        else if (strncmp(snameshort2[pj],"det",3)==0)
          fprintf(sp,", minimum=_ZERO_, no_river_dilution=detritus_no_river_dilution "); 
        else if (strncmp(snameshort2[pj],"nut",3)==0)
          fprintf(sp,", minimum=_ZERO_, no_river_dilution=nutrient_no_river_dilution "); 
        else
	 fprintf(sp,", minimum=_ZERO_, no_river_dilution=.%s. ",yesno[ss]); 
	}
      else
	if(strstr(pcom[ni][pj],"size")!=NULL)
	 fprintf(sp,", minimum=-2.0d0, no_river_dilution=.%s. ",yesno[ss]); 
	else
	 fprintf(sp,", minimum=_ZERO_, no_river_dilution=.%s. ",yesno[ss]); 
	  
      if(setvel[ni][pj]>=0)
        fprintf(sp,", vertical_movement=%s/secs_pr_day",parname[setvel_n[ni][pj]][setvel[ni][pj]]); 
//      if(trait[pj] || strstr(snameshort2[pj],"phy"))	
//      else  fprintf(sp,"%s\tvertical_movement=_ZERO_)\n",tabs[(pjs<0)],'%'); 
      fprintf(sp,")\n");
      if(SED) // fabm sediment models distinguish between particulate and dissolved states
        {
	if(strstr(pcom[ni][pj],"det")!=NULL || strstr(snameshort2[pj],"mpb")!=NULL) ss=1; else ss=0; 
        fprintf(sp,"%scall self%cset_variable_property(self%cid_%s,'particulate',.%s.)\n",indent0,'%','%',snameshort2[pj],yesno[ss]);
	}

      ss=-1;
      velem= snameshort2[pj][strlen(snameshort2[pj])-1];
      for(pi=0;pi<nelements &&ss<0;pi++)
        if(elements[pi]==velem) ss=pi;
      ChemSpec[pj]=ss;
      if(out)printf(" reg_state(self%cid_%s,%s'%s',... %c %d\n",'%',snameshort2[pj],tabs[(pjs<0)],snameshort2[pj],velem,ss);
      if(ss>=0)
       fprintf(sp,"%s%scall self%cadd_to_aggregate_variable(standard_variables%ctotal_%s,self%cid_%s)\n",indent0,tabs[(pjs<0)],'%','%',longelemnam[ss],'%',snameshort2[pj]);
      }
      
     if( pj==nump[ni]-1 && sws==1) fprintf(sp,"end if\n");
     }
if (bdim>1)
  fprintf(sp,"end do\n");
fprintf(sp,"\n!!------- Register diagnostic variables  ------- \n");
//if(strcmp(modname,"maecs")==0)  fprintf(sp,"if (self%cDebugDiagOn) then\n",'%');
ni=nir+1;
for(pjs=-1;pjs<nump[nis];pjs++)
 {
 for(pj=0,sws=0;pj<nump[ni];pj++)
   if(swi[ni][pj]==pjs)
     {
     sl[pj]= ' '; if(strlen(snameshort2[pj])<7) sl[pj]='\t';
     if(sws==0 && pjs>=0)
       fprintf(sp,"\nif (self%c%s) then\n",'%',swin[ni][pj]), sws=1;

     fprintf(sp,"%scall self%cregister_diagnostic_variable(self%cid_%s,%s'%s',",indent0,'%','%',
       parname[ni][pj],fil(sc,parname[ni][pj],8),parname[ni][pj]);
     fprintf(sp,"'%s', '%s %s', &\n%s  ",punit[ni][pj],pcom[ni][pj],parname[ni][pj],indent0);
//  fprintf(sp,"time_treatment=time_treatment_%s)\n",partypen[ni][pj]);  
//  fprintf(sp,"output=output_instantaneous)\n");  
     fprintf(sp,"output=output_time_step_averaged)\n");  
     if(out)printf(" reg_diag(self%cid_%s,'%s',...\t%s\n",'%',parname[ni][pj],parname[ni][pj],swin[ni][pj]);
     }
 if(sws==1) fprintf(sp,"end if\n");
 }  
//if(strcmp(modname,"maecs")==0)  fprintf(sp,"end if\n");
 /*
if(nelements>0)
  {
  fprintf(sp,"\n!!------- Register conserved quantities  ------- \n");
  for(pi=0;pi<nelements;pi++)
   fprintf(sp,"%scall self%cregister_conserved_quantity(self%cid_tot%c,'%c','mmol-%c/m**3','total-%c')\n",
     indent0,'%','%',elements[pi],elements[pi],elements[pi],elements[pi]);
  }
    */
pi=0;
if(strcmp(modname,"maecs")==0)
  {
  if(swi1[pi][0]!='-')  fprintf(sp,"\n! ------ check dependencies in diag switches -------\n");
  while(swi1[pi][0]!='-')
   {
   fprintf(sp,"if (self%c%s .and. .not. self%c%s) call self%cfatal_error('maecs_init','%s=TRUE requires %s=TRUE')\n",'%',swi1[pi],'%',swi2[pi],'%',swi1[pi],swi2[pi]);
   pi++;  
   }
  }
ni=nir;
fprintf(sp,"\n!!------- Register environmental dependencies  ------- \n");

for(pjs=-1;pjs<nump[nis];pjs++)
  for(pj=0,sws=0;pj<nump[ni];pj++)
     { 
     if(swi[ni][pj]==pjs)
       {      
       if(sws==0 && pjs>=0)
         fprintf(sp,"\nif (self%c%s) then\n",'%',swin[ni][pj]), sws=1;
//   fprintf(sp,"%scall self%cregister_dependency(self%cid_%s,varname_%s%s)\n",'%','%',
//       if(partypen[ni][pj][0]=='h' && strstr(parname[ni][pj],"tot")!=NULL )
       strcpy(line2,FabmDepVarName),strcpy(line,partypen[ni][pj]);  
       strcpy(line1,pmapstring[ni][pj]);
       if(strstr(parname[ni][pj],"vert")!=NULL || strstr(parname[ni][pj],"flux")!=NULL || strstr(parname[ni][pj],"_dep")!=NULL )
         {
	 strcpy(line2,"");
         if(strstr(parname[ni][pj],"diag")==NULL) 
	   strcpy(line,"dependency");
	 else
	   strcat(line1,", output=output_time_step_averaged");
	 }

       if(out) printf(" reg dep(self%cid_%s,varname_%s) \t switch=%d %d %d\n",'%',parname[ni][pj],pmapstring[ni][pj] ,swi[ni][pj],sws,pjs);
       fprintf(sp,"%s%scall self%cregister_%s(self%cid_%s,%s%s)\n",indent0,tabs[(pjs<0)],'%',line,'%',parname[ni][pj],line2,line1);
         
//    }   else     {printf("\n** ERROR: external forcing %s not found !!!\n now exit...\n",'%',parname[ni][pj]);exit(0);}
       } // if(swi[ni][pj]
     if( pj==nump[ni]-1 && sws==1) fprintf(sp,"end if\n");
     }  
     
if(strlen(init_incl)>1)
  {
  printf("\n! extra lines included from %s \n",init_incl);
  fprintf(sp,"\n! extra lines included from %s \n",init_incl);
  sp4=fopen(init_incl,"r"); 
  if(sp4==NULL) {printf("Error while opening %s !\n\n",init_incl),exit(0);}

  lr=(char *)1;
  while(lr!=NULL)
    {
    lr=fgets(line,256,sp4); 
    fputs(line,sp);    
    }
  fclose(sp4);
  }
fprintf(sp,"\n\n%sreturn\n\n!!-------  if files are not found ...  \n",indent0);
for(ni=0;ni<nir;ni++)
  fprintf(sp,"%d call self%cfatal_error('%s_init','Error reading namelist %s.')\n",90+ni,'%',modname,nmlname[ni]);
for(ni=0;ni<nir;ni++)
  fprintf(sp,"%d call self%cfatal_error('%s_init','Namelist %s was not found in file.')\n",90+ni+MI,'%',modname,nmlname[ni]);

fprintf(sp,"\nend subroutine initialize\n\n");
//fprintf(sp,"\nend module %s_initvar\n",modname);
fprintf(sp,"!!----------------------------------------------------------------------\n");
fprintf(sp,"!!   end of section generated by parser \n");
fprintf(sp,"!!----------------------------------------------------------------------\n");
//fclose(sp);

if(TYPES)
  {
//  printf("closing %ld %ld %ld %ld  ...\n",sp,sp1,sp2,sp3);  
  cp1=strstr(modname,"#SP#");lr=(char *)1;
  while(lr!=NULL && cp1==NULL)
    {
    lr=fgets(line,256,sp1); 
    cp1=strstr(line,"#SP#");
    }
  while(lr!=NULL )
    {
    fputs(line,sp);
    lr=fgets(line,256,sp1); 
    }
  cp1=strstr(modname,"#SP#");
  lr=(char *)1;
  while(lr!=NULL && cp1==NULL)
    {
    lr=fgets(line,256,sp2); 
    cp1=strstr(line,"#SP#");
    }
  while(lr!=NULL )
    {  
    fputs(line,sp3);
    lr=fgets(line,256,sp2); 
    }
    
  fclose(sp);
  fclose(sp1);
  fclose(sp2);
  fclose(sp3);

// read existing code before listed initialisation  
  sprintf(insname,"%s%s_do.F90\0",dirn_f90,modname);
  printf("reading from %s ...\n",insname);
  sp1=fopen(insname,"r"); 
  if(sp1==NULL) {printf("Problem while opening %s !\n",insname),exit(0);}
  sprintf(outn,"%s_do_gen.F90\0",modname);
  sp=fopen(outn,"w"); 
  if(sp==NULL) {printf("Problem while opening %s !\n",outn),exit(0);}
  if(out || 1) printf("writing lists into %s...\n",outn); 
  lr=fgets(line,256,sp1); 

  }
else
  {
  lr=(char *)1;
  //  printf("lr=%c\n",lr);
  cp1=strstr(modname,"#SP#");
  while(lr!=NULL && cp1==NULL)
    {
    lr=fgets(line,256,sp1); 
//    fputs(line,sp3);
    cp1=strstr(line,"#SP#");
    }
  }
  
  //sprintf(outn,"%s_gen.F90\0",modname);
//printf("re-opening %s ...\n",outn);
//sp3=fopen(outn,"a");spv[0]=sp3;
spv[0]=sp;

if(NewModF90)
  {
  sprintf(outn,"%s_ins.F90\0",modname);
  sp=fopen(outn,"w"); spv[1]=sp;
  }

/* -----------------------------------------------------------*/
/*  insert formatted lists at tag position in existing code   */
/* -----------------------------------------------------------*/
for(d=0;d<1+0*NewModF90;d++)
  {

  while(lr!=NULL )
    {
    fputs(line,spv[d]);
//  strncpy(line2,line,12);
//  printf("%s\t%d\n",line2,strlen(line));
    lr=fgets(line,256,sp1); 
    cp1=strstr(line,"!#");
    if(cp1 != NULL)
      {
      fputs(line,spv[d]);
      pc=-1;
      for(pi=0;pi<3;pi++) if (strncmp(cp1+4,keys[pi],3)==0) pc=pi;
      if(strncmp(cp1+4,"GED",3)==0) pc=5;
      if(strncmp(cp1+4,"DIA",3)==0) pc=6;   
      if(strncmp(cp1+4,"CON",3)==0) pc=7;  
      
      if(pc==6) 
        printf("\n diags: ni=%d pj<%d\n\n", nir+1,nump[ni]);  
      
      
      if(pc>=0 && pc<5)
        {
        fprintf(spv[d],"!---------- %s for each state variable ----------\n",keys[pc]);
        printf("!---------- %s for each state variable ----------\n",keys[pc]);
        for(ni=ni0,pjs=-1;pjs<nump[nis];pjs++)
         for(pj=0,sws=0;pj<nump[ni];pj++)
          {
          if(swi[ni][pj]==pjs)
           {
           if(sws==0 && pjs>=0)
             fprintf(spv[d],"if (self%c%s) then\n",'%',swin[ni][pj]), sws=1;  
//           if(trait[pj] && strlen(traitpre)>0)   strcpy(outn,traitpre);
//           else   strcpy(outn,"\0");
//           strcat(outn,snameshort[pj]);//snameshort3
	   switch(pc) {
	     case 0:
	      fprintf(spv[d],"! --- dynamics of %s ~ %s\n",snameshort[pj],pcom[ni][pj]); 
              fprintf(spv[d],"!%s rhsv%c%s = . %s%c%s\n",tabs[(pjs<0)],'%',snameshort2[pj],vstructn,'%',snameshort2[pj]);
	      break;
	     case 1:
	      if(SCALEFAC==0)  sprintf(outn," UNIT");
	      else sprintf(outn,"");
              fprintf(spv[d],"%s  _SET_ODE_(self%cid_%s, rhsv%c%s%s)\n",
//		       tabs[(pjs<0)],'%',snameshort2[pj],'%',snameshort2[pj],outn);
		       tabs[(pjs<0)],'%',idname[pj],'%',snameshort2[pj],outn);
	      break;
	     case 2:
//              fprintf(spv[d],"%s  _GET_(self%cid_%s, var%c%s)  ! %s in %s\n",tabs[(pjs<0)],
// all biomass=non-trait variables are assumed to be in special structures, if not part of
//    overall variable container

//printf("%s \t%d %d\n",snameshort2[pj],trait[pj],strcmp(modname,"maecs"));

              if(trait[pj]==0 && (strchr(vstructc,"A")==NULL || strchr(vstructc,"E")==NULL) && strcmp(modname,"maecs")==0 && ChemSpec[pj]>=0)/*strlen(vstructn)<1*/
	        {
		line2[0]='\0'; 
		strncat(line2,snameshort2[pj],strlen(snameshort2[pj])-1);
	        velem = snameshort2[pj][strlen(snameshort2[pj])-1];
	        if(velem=='S') 
                  sprintf(line,"%s%c%ci",line2,'%',velem);
	        else
                  sprintf(line,"%s%c%c",line2,'%',velem);
		}
	      else if(trait[pj]==1) // trait structure in MAECS
                sprintf(line,"phy%c%s",'%',snameshort2[pj]);
	      else
                sprintf(line,"%s%c%s",vstructn,'%',snameshort2[pj] );
		
              fprintf(spv[d],"%s  _GET_(self%cid_%s, %s)  ! %s in %s\n",tabs[(pjs<0)],
//		       '%',snameshort2[pj],line,pcom[ni][pj],punit[ni][pj]);
		       '%',idname[pj],line,pcom[ni][pj],punit[ni][pj]);
	     }
	   }
          if( pj==nump[ni]-1 && sws==1) fprintf(spv[d],"end if\n");
	  }	    
	 }
      if(pc==5) 
        for(ni=nir,pj=0;pj<nump[ni];pj++)
	  {
	  if(strchr(vstructc,'E')!=NULL || strchr(vstructc,'A')!=NULL) 
             sprintf(line,"%s%c%s",vstructn,'%',parname[ni][pj]);
	  else
             sprintf(line,"%s",parname[ni][pj]);

	  if(partypen[ni][pj][0]=='h')
            fprintf(spv[d],"  _GET_HORIZONTAL_(self%cid_%s,%s) ! %s\n",'%',parname[ni][pj],line,pcom[ni][pj]);
	  else
            fprintf(spv[d],"  _GET_(self%cid_%s, %s)  ! %s\n",'%',parname[ni][pj],line,pcom[ni][pj]);
	  }
//            fprintf(spv[d],"  _GET_(self%cid_%s, var%c%s)  ! %s\n",'%',parname[ni][pj],'%',parname[ni][pj],pcom[ni][pj]);
      if(pc==6) 
	for(ni=nir+1,pjs=-1;pjs<nump[nis];pjs++)
	 { 
         for(pj=0,sws=0;pj<nump[ni];pj++)
	  if(strstr(parname[ni][pj],"tot")==NULL ) // total nut assumed horizontal
           if(swi[ni][pj]==pjs)
            {
            if(sws==0 && pjs>=0)
              fprintf(spv[d],"if (self%c%s) then\n",'%',swin[ni][pj]), sws=1;  
	    sprintf(line2,"%s, _REPLNAN_(%s",parname[ni][pj], pmapstring[ni][pj]);
	    if(found_rate[ni][pj]==1 && SCALEFAC && strstr(pcom[ni][pj],"RHS")!=NULL )  strcat(line2,"*secs_pr_day");  
            fprintf(spv[d],"  _SET_DIAGNOSTIC_(self%cid_%s",'%',line2);
            fprintf(spv[d],"))%s!%s %s\n",fil(sc,line2,32),partypen[ni][pj],pcom[ni][pj]);
	    }
	  if(sws==1) fprintf(sp,"end if\n");
	  } 
      if(pc==7) /* get_conserved_quantities */
        {
        fprintf(spv[d]," real(rk) :: ");
        for(ni=ni0,pj=0;pj<nump[ni];pj++)
           fprintf(spv[d]," %s%c",snameshort2[pj],tcb[pj==nump[ni]-1]); 

	for(ni=ni0,pj=0;pj<nump[ni];pj++)
           fprintf(spv[d]," _GET_(self%cid_%s, %s)  ! %s in %s\n",
	    '%',snameshort2[pj],snameshort2[pj],pcom[ni][pj],punit[ni][pj]);
	  
         printf("nelements=%d  \n",nelements);
/*       for(pj=0;pj<nelements;pj++)
            fprintf(spv[d]," _SET_CONSERVED_QUANTITY_(self%cid_tot%c,   )\n",'%',elements[pj]);
        for(pj=0;pj<nelements;pj++)
        printf("%c  \n",elements[pj]);
	  */
	  }
      cp1=strstr(modname,"#SP#");
      while(lr!=NULL  && cp1==NULL)
        {
        lr=fgets(line,256,sp1); 
        printf("delete %s",line);
        cp1=strstr(line,"!#E");
        } 	 
      } 
    }

  if (out) printf("closing file %d...\n",d);
  fclose(spv[d]);
  }
//fclose(sp);
//fclose(sp3);
fclose(sp1);

return;
}

char *fil(char *sc,char *pre, int maxn)
{
int pl;
char spc[49]="                                        ";

//printf("%s...%d\n",pre,strlen(pre));
pl=maxn-strlen(pre);
if(pl<1) pl=1;
strcpy(sc,spc);
sc[pl]='\0';sc[pl+1]='\0';
//printf("%s...\n",sc);
return sc;
}

char *insul(char *str)
{
char *pch,spc[99]="";

pch=strchr(str,'_');
while (pch!=NULL)
  {
  strncat(spc,str,pch-str);
  strcat(spc,"\\_");
  str=pch+1;
  pch=strchr(str,'_');
  }
  
strcat(spc,str);
//printf("%s...\n",sc);
return spc;
}




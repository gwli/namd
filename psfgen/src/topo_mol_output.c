
#include <stdlib.h>
#include "topo_mol_output.h"
#include "topo_mol_struct.h"
#include "pdb_file.h"

int topo_mol_write_pdb(topo_mol *mol, FILE *file, void *v, 
                                void (*print_msg)(void *, const char *)) {

  int iseg,nseg,ires,nres,atomid;
  int has_guessed_atoms = 0;
  double x,y,z,o,b;
  topo_mol_segment_t *seg;
  topo_mol_residue_t *res;
  topo_mol_atom_t *atom;

  if ( ! mol ) return -1;

  write_pdb_remark(file,"original generated coordinate pdb file");

  atomid = 0;
  nseg = hasharray_count(mol->segment_hash);
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);
      for ( atom = res->atoms; atom; atom = atom->next ) {
        /* Paranoid: make sure x,y,z,o are set. */
        x = y = z = 0.0; o = -1.0;
        ++atomid;
        switch ( atom->xyz_state ) {
        case TOPO_MOL_XYZ_SET:
          x = atom->x;  y = atom->y;  z = atom->z;  o = 1.0;
          break;
        case TOPO_MOL_XYZ_GUESS:
        case TOPO_MOL_XYZ_BADGUESS:
          x = atom->x;  y = atom->y;  z = atom->z;  o = 0.0;
          has_guessed_atoms = 1;
          break;
        default:
          print_msg(v,"ERROR: Internal error, atom has invalid state.");
          print_msg(v,"ERROR: Treating as void.");
          /* Yes, fall through */
        case TOPO_MOL_XYZ_VOID:
          x = y = z = 0.0;  o = -1.0;
          break;
        }
        b = atom->partition;
        write_pdb_atom(file,atomid,atom->name,res->name,atoi(res->resid),
		"",(float)x,(float)y,(float)z,(float)o,(float)b,"",seg->segid,
		atom->element);
      }
    }
  }

  write_pdb_end(file);
  if (has_guessed_atoms) {
    print_msg(v, 
        "Info: Atoms with guessed coordinates will have occupancy of 0.0.");
  }
  return 0;
}

int topo_mol_write_psf(topo_mol *mol, FILE *file, int charmmfmt, void *v, 
                                void (*print_msg)(void *, const char *)) {

  char buf[128];
  int iseg,nseg,ires,nres,atomid;
  topo_mol_segment_t *seg;
  topo_mol_residue_t *res;
  topo_mol_atom_t *atom;
  topo_mol_bond_t *bond;
  int nbonds;
  topo_mol_angle_t *angl;
  int nangls;
  topo_mol_dihedral_t *dihe;
  int ndihes;
  topo_mol_improper_t *impr;
  int nimprs;
  int numinline;

  if ( ! mol ) return -1;

  fprintf(file,"PSF\n\n%8d !NTITLE\n",1);
  if ( charmmfmt ) 
   fprintf(file," REMARKS %s\n","original generated structure charmm psf file");
  else
   fprintf(file," REMARKS %s\n","original generated structure x-plor psf file");
  fprintf(file,"\n");

  atomid = 0;
  nbonds = 0;
  nangls = 0;
  ndihes = 0;
  nimprs = 0;
  nseg = hasharray_count(mol->segment_hash);
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);
      for ( atom = res->atoms; atom; atom = atom->next ) {
        atom->atomid = ++atomid;
        for ( bond = atom->bonds; bond;
                bond = topo_mol_bond_next(bond,atom) ) {
          if ( bond->atom[0] == atom && ! bond->del ) {
            ++nbonds;
          }
        }
        for ( angl = atom->angles; angl;
                angl = topo_mol_angle_next(angl,atom) ) {
          if ( angl->atom[0] == atom && ! angl->del ) {
            ++nangls;
          }
        }
        for ( dihe = atom->dihedrals; dihe;
                dihe = topo_mol_dihedral_next(dihe,atom) ) {
          if ( dihe->atom[0] == atom && ! dihe->del ) {
            ++ndihes;
          }
        }
        for ( impr = atom->impropers; impr;
                impr = topo_mol_improper_next(impr,atom) ) {
          if ( impr->atom[0] == atom && ! impr->del ) {
            ++nimprs;
          }
        }
      }
    }
  }
  sprintf(buf,"total of %d atoms",atomid);
  print_msg(v,buf);
  sprintf(buf,"total of %d bonds",nbonds);
  print_msg(v,buf);
  sprintf(buf,"total of %d angles",nangls);
  print_msg(v,buf);
  sprintf(buf,"total of %d dihedrals",ndihes);
  print_msg(v,buf);
  sprintf(buf,"total of %d impropers",nimprs);
  print_msg(v,buf);

  fprintf(file,"%8d !NATOM\n",atomid);
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);
      if ( charmmfmt ) for ( atom = res->atoms; atom; atom = atom->next ) {
        fprintf(file,"%8d %-4s %-4s %-4s %-4s %4d %10.6f     %9.4f  %10d\n",
                atom->atomid, seg->segid,res->resid,res->name,
                atom->name,atom->typeid,atom->charge,atom->mass,0);
      } else for ( atom = res->atoms; atom; atom = atom->next ) {
        fprintf(file,"%8d %-4s %-4s %-4s %-4s %-4s %10.6f     %9.4f  %10d\n",
                atom->atomid, seg->segid,res->resid,res->name,
                atom->name,atom->type,atom->charge,atom->mass,0);
      }
    }
  }
  fprintf(file,"\n");

  fprintf(file,"%8d !NBOND: bonds\n",nbonds);
  numinline = 0;
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);
      for ( atom = res->atoms; atom; atom = atom->next ) {
        for ( bond = atom->bonds; bond;
                bond = topo_mol_bond_next(bond,atom) ) {
          if ( bond->atom[0] == atom && ! bond->del ) {
            if ( numinline == 4 ) { fprintf(file,"\n");  numinline = 0; }
            fprintf(file," %7d %7d",atom->atomid,bond->atom[1]->atomid);
            ++numinline;
          }
        }
      }
    }
  }
  fprintf(file,"\n\n");

  fprintf(file,"%8d !NTHETA: angles\n",nangls);
  numinline = 0;
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);
      for ( atom = res->atoms; atom; atom = atom->next ) {
        for ( angl = atom->angles; angl;
                angl = topo_mol_angle_next(angl,atom) ) {
          if ( angl->atom[0] == atom && ! angl->del ) {
            if ( numinline == 3 ) { fprintf(file,"\n");  numinline = 0; }
            fprintf(file," %7d %7d %7d",atom->atomid,
                angl->atom[1]->atomid,angl->atom[2]->atomid);
            ++numinline;
          }
        }
      }
    }
  }
  fprintf(file,"\n\n");

  fprintf(file,"%8d !NPHI: dihedrals\n",ndihes);
  numinline = 0;
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);
      for ( atom = res->atoms; atom; atom = atom->next ) {
        for ( dihe = atom->dihedrals; dihe;
                dihe = topo_mol_dihedral_next(dihe,atom) ) {
          if ( dihe->atom[0] == atom && ! dihe->del ) {
            if ( numinline == 2 ) { fprintf(file,"\n");  numinline = 0; }
            fprintf(file," %7d %7d %7d %7d",atom->atomid,
                dihe->atom[1]->atomid,dihe->atom[2]->atomid,
                dihe->atom[3]->atomid);
            ++numinline;
          }
        }
      }
    }
  }
  fprintf(file,"\n\n");

  fprintf(file,"%8d !NIMPHI: impropers\n",nimprs);
  numinline = 0;
  for ( iseg=0; iseg<nseg; ++iseg ) {
    seg = mol->segment_array[iseg];
    if (! seg) continue;
    nres = hasharray_count(seg->residue_hash);
    for ( ires=0; ires<nres; ++ires ) {
      res = &(seg->residue_array[ires]);
      for ( atom = res->atoms; atom; atom = atom->next ) {
        for ( impr = atom->impropers; impr;
                impr = topo_mol_improper_next(impr,atom) ) {
          if ( impr->atom[0] == atom && ! impr->del ) {
            if ( numinline == 2 ) { fprintf(file,"\n");  numinline = 0; }
            fprintf(file," %7d %7d %7d %7d",atom->atomid,
                impr->atom[1]->atomid,impr->atom[2]->atomid,
                impr->atom[3]->atomid);
            ++numinline;
          }
        }
      }
    }
  }
  fprintf(file,"\n\n");

  fprintf(file,"%8d !NDON: donors\n",0);
  fprintf(file,"\n\n");

  fprintf(file,"%8d !NACC: acceptors\n",0);
  fprintf(file,"\n\n");

  fprintf(file,"%8d !NNB\n\n",0);
  /* Pad with zeros, one for every atom */
  {
    int i, fullrows;
    fullrows = atomid/8;
    for (i=0; i<fullrows; ++i) 
      fprintf(file, "%8d%8d%8d%8d%8d%8d%8d%8d\n",0,0,0,0,0,0,0,0);
    for (i=atomid - fullrows*8; i; --i)
      fprintf(file, "%8d",0);
  } 
  fprintf(file,"\n\n");

  fprintf(file,"%8d %7d !NGRP\n%8d%8d%8d\n",1,0,0,0,0);
  fprintf(file,"\n");

  return 0;
}


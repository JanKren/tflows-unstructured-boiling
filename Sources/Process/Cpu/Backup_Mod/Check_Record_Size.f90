!==============================================================================!
  subroutine Check_Record_Size(Backup, Comm, var_name, vs_file, vs_here)
!------------------------------------------------------------------------------!
!   Compares the record size stored in the backup file against the size this   !
!   grid would write for the same record, and stops the run when they differ.  !
!                                                                              !
!   Every Save_* routine writes the payload size in bytes right after the      !
!   record name, and computes that size as a pure function of the grid and     !
!   the working precision - (nc_tot + nb_tot) * RP for a cell field,           !
!   (3*nc_tot + 3*nb_tot) * RP for a variable, length * IP for an array, and   !
!   so on.  Until now every Load_* routine read that word back and used it     !
!   only to skip records it did not want, never to check the record it was     !
!   about to read.  A backup written on one grid was therefore read into a     !
!   different grid's arrays without a word of complaint.                       !
!                                                                              !
!   The case that motivated this is a restart after grid adaptation.           !
!   Rebuild_After_Change calls Form_Maps_For_Backup after every sweep, so the  !
!   backup is written in refined-grid global numbering, while the .cfn and     !
!   .dim on disk still describe the base grid and nothing ever rewrites them.  !
!   The restart then loads the base grid and reads the refined-grid image      !
!   into it.  The record search still succeeds, because it matches on the      !
!   name; the payload read then uses the current grid's nc_tot.  The result    !
!   is a plausible-looking field made of misaligned bytes.                     !
!                                                                              !
!   The check costs one integer comparison per record actually read, changes   !
!   nothing on disk, and works on every backup T-Flows has ever written,       !
!   because the size word has always been there.  Besides adaptation it also   !
!   catches a restart against a regenerated or re-divided mesh, and a          !
!   REAL=single binary pointed at a double-precision backup.                   !
!                                                                              !
!   Collectiveness.  Comm % Read_Int sets a plain shared view and every rank   !
!   reads the same four bytes, so vs_file is bit-identical everywhere and all  !
!   ranks take the same branch.  Message % Error calls Global % End_Parallel   !
!   and stops, which is safe only under that property.  Keep it that way: a    !
!   record size that could differ between ranks would turn this stop into a    !
!   hang.                                                                      !
!------------------------------------------------------------------------------!
  implicit none
!---------------------------------[Arguments]----------------------------------!
  class(Backup_Type) :: Backup
  type(Comm_Type)    :: Comm
  character(len=*)   :: var_name
  integer            :: vs_file  ! record size as stored in the file
  integer            :: vs_here  ! record size this grid would write
!-----------------------------------[Locals]-----------------------------------!
  character(1024) :: text
!------------------------[Avoid unused parent warning]-------------------------!
  Unused(Backup)
!==============================================================================!

  if(vs_file .eq. vs_here) return

  write(text, '(a,a,a,   i0,a,i0,a,   i0,a,i0,a)')                          &
    'Backup record "', trim(var_name), '" does not fit this grid.  The '//  &
    'file holds ', vs_file, ' bytes for it, this grid would write ',        &
    vs_here, ' bytes.  This grid has ',                                     &
    Comm % nc_tot, ' interior and ', Comm % nb_tot, ' boundary cells.  '//  &
    'The backup was written on a different grid, so reading it would '//    &
    'silently fill the fields with misaligned values.  Exiting!'

  call Message % Error(72, trim(text), file=__FILE__, line=__LINE__)

  end subroutine

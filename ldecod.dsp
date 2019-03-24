# Microsoft Developer Studio Project File - Name="ldecod" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

CFG=ldecod - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "ldecod.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "ldecod.mak" CFG="ldecod - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "ldecod - Win32 Release" (based on "Win32 (x86) Console Application")
!MESSAGE "ldecod - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "ldecod - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "ldecod\Release"
# PROP Intermediate_Dir "ldecod\Release"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /Ob2 /I "ldecod\inc" /I "lcommon\inc" /D "WIN32" /D "NDEBUG" /D "_CONSOLE" /D "_MBCS" /Fr /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /machine:I386 /out:"./bin/ldecod.exe"

!ELSEIF  "$(CFG)" == "ldecod - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "ldecod___Win32_Debug"
# PROP BASE Intermediate_Dir "ldecod___Win32_Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "ldecod/Debug"
# PROP Intermediate_Dir "ldecod/Debug"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /GX /Zi /Od /I "ldecod/inc" /I "lcommon/inc" /D "WIN32" /D "_DEBUG" /D "_CONSOLE" /D "_MBCS" /FR /YX /FD /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /debug /machine:I386 /pdbtype:sept
# ADD LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /subsystem:console /profile /debug /machine:I386 /out:"./bin/ldecod.exe"

!ENDIF 

# Begin Target

# Name "ldecod - Win32 Release"
# Name "ldecod - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Source File

SOURCE=.\ldecod\src\abt.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\abt_const.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\b_frame.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\biaridecod.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\bitsbuf.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\block.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\cabac.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\decodeiff.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\erc_api.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\erc_do_i.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\erc_do_p.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\errorconcealment.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\filehandle.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\fmo.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\golomb_dec.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\header.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\image.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\ldecod.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\leaky_bucket.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\loopFilter.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\macroblock.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\mbuffer.c
# End Source File
# Begin Source File

SOURCE=.\lcommon\src\memalloc.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\nal.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\nal_part.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\output.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\rtp.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\sei.c
# End Source File
# Begin Source File

SOURCE=.\ldecod\src\uvlc.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=.\ldecod\inc\abt.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\b_frame.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\biaridecod.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\bitsbuf.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\block.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\cabac.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\contributors.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\decodeiff.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\defines.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\elements.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\erc_api.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\erc_do.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\erc_globals.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\errorconcealment.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\fmo.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\global.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\golomb_dec.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\header.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\image.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\leaky_bucket.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\macroblock.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\mbuffer.h
# End Source File
# Begin Source File

SOURCE=.\lcommon\inc\memalloc.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\rtp.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\sei.h
# End Source File
# Begin Source File

SOURCE=.\ldecod\inc\uvlc.h
# End Source File
# End Group
# Begin Source File

SOURCE=.\bin\decoder.cfg
# End Source File
# End Target
# End Project

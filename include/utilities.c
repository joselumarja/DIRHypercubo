#include "utilities.h"
#include "definitions.h"

#include <stdio.h>
#include <errno.h>
#include <string.h>
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>

int isFile(const char *path)
{
	FILE *FilePointer;
	if((FilePointer = fopen(path,"r"))==NULL)
	{
		fprintf(stderr, "Error in open file %s : %s\n",path,strerror(errno));
		return false;
	}
	
	fclose(FilePointer);
	
	return true;
}

int countDataNumber(const char *FilePath)
{
	FILE *FilePointer;
	
	int Count=0,i;
	
	char Line[MAX_DATA_BLOCK];
	
	if((FilePointer = fopen(FilePath,"r"))==NULL)
	{
		fprintf(stderr, "Error in file %s : %s\n",FilePath,strerror(errno));
		return false;
	}
	
	while(fgets(Line,MAX_DATA_BLOCK-1,FilePointer)!=NULL)
	{
		for(i=0; i<MAX_DATA_BLOCK || Line[i]=='\0'; i++)
		{
			if(Line[i]==',')
			{
				Count++;
			}
		}
	}
	
	fclose(FilePointer);
	
	return Count;
}
		
			

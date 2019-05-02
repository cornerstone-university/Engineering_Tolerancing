/*
Name: Tolerance Analysis Program
Author: Avery Wittmer
Desc.: This program parses an inputted tolerance table and outputs a Basic Gap analysis.
It then builds on its analysis through suggesting ways to close the gap and dealing with random dimensions.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


void Calculation();
void Method_2(int, float*, int*, float*, char*, float*, float*, float, float);
void Method_3(int, float*, int*, float*, char*, float, float);
void random_dimension(float, float, float *);
void MonteCarlo(int, float*, int*, float*, char*, float, float*, float*);
float Mean(float*, float*);
float Deviation(float*, float, float*);


int main()
{
    int s = 100, n = 0, sign[s], numParts = 0, numGaps = 0, numLines = 0, scan = 0, input, aboveMax = 0, belowMin = 0;
    FILE* fpin;
    FILE* fpout;
    FILE* csv;
    char fv[s], file[s];
    float size[s], tol[s], minGap[s], maxGap[s], totalMinG = 0.0, totalMaxG = 0.0, absMinGap, absMaxGap,
    correctGap, correctTol, bilatGap = 0.0, bilatTol = 0.0, randomSize, randomGap, randomGaps[10000],
    mean, dev;

if((fpin=fopen("ToleranceTable.txt", "r")) == NULL){
    printf("File DNE.\n");
    exit(-1);
}

fpout = fopen("Output.txt", "w");
csv = fopen("ExcelOutput.csv", "w");

if(fgets(file, sizeof(file), fpin) == NULL){
    printf("File is empty.\n");
    exit(-1);
}
rewind(fpin);

for(int i=0; fscanf(fpin, "%c", &file[n]) != EOF; i++){
    n++;
}
rewind(fpin);

for(int i=0; i<n; i++){
    if(file[i] == 'A'){
        if(file[i-1] == 'P'){

        scan = fscanf(fpin, "PART,%f,%d,%f,%c\n", &size[numParts], &sign[numParts], &tol[numParts], &fv[numParts]);
        fprintf(csv, "%.3f, ", size[numParts]);
        numParts++;
    }
    else if(file[i-1] == 'G'){

        scan = fscanf(fpin, "GAP,%f,%f\n", &minGap[numGaps], &maxGap[numGaps]);
        totalMinG += minGap[numGaps];
        totalMaxG += maxGap[numGaps];
        numGaps++;
        }
    }
}
rewind(fpin);

absMaxGap = maxGap[0];
absMinGap = minGap[0];
for(int i=0; i<numGaps; i++){
    if(minGap[i] < absMinGap){
        absMinGap = minGap[i];
    }
    else if(maxGap[i] > absMaxGap){
        absMaxGap = maxGap[i];
    }
}

numLines = numParts + numGaps;
correctGap = (totalMinG + totalMaxG) / 2;
correctTol = (totalMinG - totalMaxG) / 2;

printf("\nWhich method of analysis do you choose? (Enter 1, 2, or 3)...\n");
printf("\nMethod 1: Base analysis of part and gap dimensions.\n");
printf("\nMethod 2: Adjust parts individually using\nnew size values for each variable part.\n");
printf("\nMethod 3: Adjust parts as a group using a tolerance scaling factor.\n\n");

scanf("%d", &input);

if(input == 1){
Calculation(numParts, numGaps, size, sign, tol, fv, minGap, maxGap, &bilatGap, &bilatTol);
}
else if(input == 2){
Method_2(numParts, size, sign, tol, fv, &bilatGap, &bilatTol, correctGap, correctTol);
}
else if(input == 3){
Method_3(numParts, size, sign, tol, fv, correctGap, correctTol);
}
else{
printf("Could not process input. Please enter 1, 2, or 3.\n");
}


printf("\nTHIRD STAGE: Monte Carlo Simulation\n");

for(int i=0; i<10000; i++){
    MonteCarlo(numParts, size, sign, tol, fv, correctGap, &randomSize, &randomGap);
    randomGaps[i] = randomGap;
    fprintf(fpout, "%.3f\n", randomGaps[i]);

    if(randomGaps[i] > absMaxGap){
        aboveMax++;
    }
    else if(randomGaps[i] < absMinGap){
        belowMin++;
    }
}
printf("\nNumber of gaps over Maximum: %d\n", aboveMax);
printf("Number of gaps under Minimum: %d\n", belowMin);

float meanVal =  Mean(randomGaps, &mean);
printf("\nMean value is %f\n", Mean(randomGaps, &mean));
printf("Standard Deviation if %f\n", Deviation(randomGaps, meanVal, &dev));


fclose(fpin);
fclose(fpout);
return 0;
}


void Calculation(int NumParts, int NumGaps, float Size[100], int Sign[100], float Tol[100], char Fv[100],
float MinGap[100], float MaxGap[100], float*BilatGap, float*BilatTol){

float MaxG = 0.0, MinG = 0.0;
int i;
FILE* CSV;
CSV = fopen("ExcelOutput.csv", "w");
     for(i=0; i<NumParts; i++){
       *BilatGap += Size[i] * Sign[i];
       *BilatTol += Tol[i];
       MaxG += Size[i] * Sign[i] + Tol[i];
       MinG += Size[i] * Sign[i] - Tol[i];
     }
       printf("Basic Gap: %.3f\"\n", *BilatGap);
       printf("Basic Gap Tolerance: %.3f\"\n", *BilatTol);
       fprintf(CSV, "%.3f, ", *BilatGap);
       fprintf(CSV, "%.3f, ", *BilatTol);

    for(i=0; i<NumGaps; i++){
    if(MaxG < MaxGap[i]){
        printf("The Maximum Gap (%.3f\") is (Less) than specified (%.3f\")\n", MaxG, MaxGap[i]);
    }
    if(MinG > MinGap[i]){
        printf("The Minimum Gap (%.3f\") is (Greater) than specified (%.3f\")\n", MinG, MinGap[i]);
    }
    if (MaxG > MaxGap[i]){
       printf("The Maximum Gap (%.3f\") is (Greater) than specified (%.3f\")\n", MaxG, MaxGap[i]);
   }
    if (MinG < MinGap[i]){
       printf("The Minimum Gap (%.3f\") is (Less) than specified (%.3f\")\n", MinG, MinGap[i]);
    }
}
}


void Method_2(int NumParts2, float Size2[100], int Sign2[100], float Tol2[100], char Fv2[100],
float *BilatGap2, float *BilatTol2, float CorrectGap, float CorrectTol){

float sizeDiff = CorrectGap - *BilatGap2;
float tolDiff = CorrectTol - *BilatTol2;
float newSize[20], newTol[20];

for(int i=0; i<NumParts2; i++){
    if(Fv2[i] == 'V'){
        newSize[i] = fabs(Size2[i] + (Sign2[i] * sizeDiff));
        newTol[i] = fabs(Tol2[i] + tolDiff);
        printf("Part #%d:\nSize = %.3f\nTolerance = %.3f\n\n", i+1, newSize[i], newTol[i]);
    }
    else{
        printf("\nPart #%d is Fixed, it cannot be changed.\n\n", i+1);
    }
  }
}


void Method_3(int NumParts3, float Size3[100], int Sign3[100], float Tol3[100], char Fv3[100],
float Correct_Gap, float Correct_Tol){

    float NewSize[20], NewTol[20], sizeVals[4], tolV=0, tolF=0, sizeFactor, tolFactor;

    for(int i=0; i<4; i++){
        sizeVals[i] = 0.0;
    }

    for(int i=0; i<NumParts3; i++){
        if(Sign3[i]>0){
            if(Fv3[i] == 'F'){
                tolF += Tol3[i];
                sizeVals[0] += Size3[i];
            }
            else if(Fv3[i] == 'V'){
                tolV += Tol3[i];
                sizeVals[1] += Size3[i];
            }
        }
        else if(Sign3[i]<0){
            if(Fv3[i] == 'F'){
                tolF += Tol3[i];
                sizeVals[2] += Size3[i];
            }
            else if(Fv3[i] == 'V'){
                tolV += Tol3[i];
                sizeVals[3] += Size3[i];
            }
        }

    }
    if(sizeVals[2] == 0){
        sizeFactor = sizeVals[3] / (sizeVals[0] - sizeVals[2] - Correct_Gap);
    }
    else{
        float a=sizeVals[1], b=sizeVals[0] - sizeVals[2] - Correct_Gap, c=sizeVals[3] * -1;

        sizeFactor=(-b+sqrtf(powf(b,2.0)-4*a*c))/(2*a);
    }

    tolFactor = (Correct_Gap - tolF) / tolV;

    for(int i=0; i<NumParts3; i++){
        if(Fv3[i] == 'V'){
            NewSize[i] = powf(sizeFactor, Sign3[i]) * Size3[i];
            NewTol[i] = tolFactor * Tol3[i];
            printf("Part #%d:\nSize = %.3f\nTolerance = %.3f\n\n", i+1, NewSize[i], NewTol[i]);
        }
        else{
           printf("\nPart #%d is Fixed, it cannot be changed.\n\n", i+1);
        }
    }
}


void random_dimension(float nominal, float tolerance, float *random_value){
	double r1, r2, r12;
	double sigma = tolerance / 3;
	do{
	r1 = (double)( rand() % 10001 ) / 10000;
	}while(r1==0);
	r2 = (double)( rand() % 10001 ) / 10000;
	r12 = sqrt(-2*log(r1))*cos(2*M_PI*r2);
	*random_value = nominal + sigma * r12;
}

void MonteCarlo(int num_Parts, float size2[100], int sign2[100], float tol2[100], char fv2[100],
float correct_Gap, float *random_Size, float *random_Gap){

    *random_Gap = 0.0;
    float randoArray[10000];

for(int i=0; i<num_Parts; i++){
    random_dimension(size2[i], tol2[i], random_Size);
    randoArray[i] = *random_Size;
    *random_Gap += randoArray[i] * sign2[i];
    }
}

float Mean(float RandomGaps[10000], float*Mean){
float sum = 0.0;
for(int i=0; i<10000; i++){
    sum += RandomGaps[i];
    }
    *Mean = sum / 10000;
    return(*Mean);
}


float Deviation(float RandGaps[10000], float Mean2, float*Dev){
float dev1 = 0.0;
for(int i=0; i<10000; i++){
    dev1 += powf(RandGaps[i] - Mean2, 2.0);
    }
    *Dev = sqrtf(dev1 / 10000);
    return(*Dev);
}

fvqe(input,codebook,codebook_size,dimension,distortion)
  float *input, **codebook;
  int codebook_size, dimension;
  float *distortion;
{
  float dist;
  int i, j, index;

  *distortion = 10000*dimension;
  index = 0;

  for (i = 0; i< codebook_size; i++)
  {
   dist = 0;
   for(j=0; j<dimension; j++)
    dist += (input[j] - codebook[i][j])*(input[j]-codebook[i][j]);
   if(dist < *distortion)
   {
    *distortion = dist;
    index = i;
   }
  }

   return index;
}
    

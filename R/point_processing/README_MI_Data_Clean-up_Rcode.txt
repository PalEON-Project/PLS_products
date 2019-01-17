MI_Data_Clean-up_Rcode.R README 
This was code used for quality control checks for both southern and northern Michigan corners. The following is a list of the steps taken with this code.
1)	Check for duplicate tree information in the northern Michigan corners
2)	Created a file with the unique FIDs (unique numerical code created by GIS) labeled for the Upper Peninsula and Lower Peninsula which Charlie used for his correction factors
3)	The recnum values continually get truncated when saving the csvs or when importing into GIS. Created code to e to add the recnum values back in using the FIDs
4)	Checekd the northern Michigan Level 0 to Level 3a taxa conversion to see what additional conversions were needed
5)	Checked that corners from townships >=11N in Tuscola, Sanilac and Huron, and >=18N in Bay, Arenac, Gladwin and Clare had the correct labels so Charlie’s correction factors are applied to the right corners
6)	labeled corners as ¼ section/section for southern Michigan 
7)	checked for duplicate southern Michigan tree information 
8)	explored southern Michigan diameter and distance values to see if they made sense.  This led to a conversation with Jack and Simon that the distances and diameters were transposed in the conversion from the mylars to the original GIS layer sent by Simon.  
9)	Checekd the southern Michigan Level 0 to Level 3a taxa conversion to see what additional conversions were needed
10)	Southern Michigan had some corners with 3-4 trees which needed to be corrected to 2 trees.  There were 68 interior corners with 4 trees which Charlie checked and corrected.  But there werew >4000 exterior townships corners with 3-4 trees.  The rule Charlie came up with for the exterior corners was: For S border corners we keep trees in quadrats 1 and 4, for E border corners we keep trees in quadrats 3 and 4.  Used R code to determine what quadrats trees are in at each corner. Labeled the corners that fell into this pattern and could be fixed easily and labeled the corners with multiple trees in the quadrats 1, 3 or 4 that means we need to re-check them in the PLS notes (e.g. on the S border a corner could have two trees in quadrat 1 and two trees in quadrat 4, so we don’t know what trees to keep without going back to the notes).

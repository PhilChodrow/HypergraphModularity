Cosponsorship Network Data

Developed by
James H. Fowler
University of California, San Diego
email: jhfowler@ucsd.edu

This file contains data on Cosponsorships in the U.S. Senate and U.S. House of Representatives for the 93rd to 108th Congresses. I have spent a lot of time developing this data so if you want to use it please cite the following two papers:

   1.     Legislative Cosponsorship Networks in the U.S. House and Senate, James H. Fowler, Social Networks 28 (4): 454-465 (October 2006)
   2.     Connecting the Congress: A Study of Cosponsorship Networks, James H. Fowler, Political Analysis 14 (4): TBD (Fall 2006)
      (Mentioned in the Washington Post, 4/11/05, and Washington Times, 4/13/05) 

These papers describe basic features of the data and how it was retrieved and processed. Because it is a very large data set there will inevitably be some mistakes in it. I do my best to continue improving the data--if you find mistakes please let me know and I will fix them.

DATA FILES

In these files an "NA" or a null value indicate the data is missing or was not matched. These occur because data was not available (e.g. early cosponsorship dates were not available at the time these files were generated) or there was a typo or other problem with the matching procedure.

    *      The SH.csv file is an aggregate results file. It is space delimited and here's a quick key:

      congress 	93rd - 108th (i.e. 1973-2004)
      chamber	House, Senate
      labels	name
      ids	ICPSR id
      ideol1	Poole and Rosenthal Common Space Score, 1st dim
      ideol2	Poole and Rosenthal Common Space Score, 2nd dim
      party	100=Dem 200=Rep
      seniority	number of Congresses served
      sponsored	number of bills sponsored
      pl	number of sponsored bills that became law
      pb	number of sponsored bills that passed chamber
      pa	number of sponsored amendments that passed chamber
      inu	unique inward cosponsors
      inb	total inward cosponsor signatures
      inw	total weighted inward cosponsor signatures (weights as described in my paper)
      outu	unique outward cosponsors
      outb	total outward cosponsor signatures
      outw	total weighted outward cosponsor signatures
      between	betweenness (unweighted)
      closeness	closeness measure -- assumes bilateral ties
      evcent	eigenvector centrality measure
      connectedness   	connectedness measure
      cc	individual clustering coefficient

The next files are all 283,994 element vectors with measures on each bill.

    *      The bills.vec file is the name of each bill as identified in the Thomas database. The name identifies the type, chamber, Congress, and number of each bill. Here's a key:

      HC   	House Concurrent Resolutions
      HE	House Resolutions
      HJ	House Joint Resolutions
      HR	House Bills
      HZ	House Amendments
      SC	Senate Concurrent Resolutions
      SE	Senate Resolutions
      SJ	Senate Joint Resolutions
      SN	Senate Bills
      SP	Senate Amendments

    *      The senate.csv file and house.csv are csv files that match ICPSR numbers ("id") to names ("name") and a few other variables for all congresses. ICPSR numbers are derived from http://voteview.com/icpsr.htm and change if a person switches party, so it is important to match by congress.

    *      The sponsors.vec file identifies the ICPSR code of each bill sponsor

    *      The cosponsors.vec file identifies the ICPSR codes of each cosponsor (one bill per line, each cosponsor is space delimited) -- large (13M)

    *      The cospcount.vec file is the total number of cosponsors on each bill

    *      The dates.vec file is the date each bill was introduced

    *      The cosponsordates.vec file shows the space delimited date(s) each bill was cosponsored -- the order of dates on each line conforms to the order of cosponsors on each line in the cosponsors.vec file -- large (22M)

    *      The party.vec file shows the party of sponsor

    *      The passedam.vec file shows whether amendment passed on the floor

    *      The passedbills.vec file shows whether bill passed on the floor

    *      The publaws.vec file shows whether bill became public law

    *      The pvtbills.vec file shows whether bill is a "private" bill 

Last Updated 1 September 2006
Copyright © 1998-2006 James Fowler, All Rights Reserved. 
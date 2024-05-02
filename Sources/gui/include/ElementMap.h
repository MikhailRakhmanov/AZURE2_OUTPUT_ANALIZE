#ifndef ELEMENTMAP_H
#define ELEMENTMAP_H

#include <QString>
#include <map>

static const std::pair<int, QString> elements[] = {
  std::pair<int, QString>(0,"n"),
  std::pair<int, QString>(1,"H"),
  std::pair<int, QString>(2,"He"),
  std::pair<int, QString>(3,"Li"),
  std::pair<int, QString>(4,"Be"),
  std::pair<int, QString>(5,"B"),
  std::pair<int, QString>(6,"C"),
  std::pair<int, QString>(7,"N"),
  std::pair<int, QString>(8,"O"),
  std::pair<int, QString>(9,"F"),
  std::pair<int, QString>(10,"Ne"),
  std::pair<int, QString>(11,"Na"),
  std::pair<int, QString>(12,"Mg"),
  std::pair<int, QString>(13,"Al"),
  std::pair<int, QString>(14,"Si"),
  std::pair<int, QString>(15,"P"),
  std::pair<int, QString>(16,"S"),
  std::pair<int, QString>(17,"Cl"),
  std::pair<int, QString>(18,"Ar"),
  std::pair<int, QString>(19,"K"),
  std::pair<int, QString>(20,"Ca"),
  std::pair<int, QString>(21,"Sc"),
  std::pair<int, QString>(22,"Ti"),
  std::pair<int, QString>(23,"V"),
  std::pair<int, QString>(24,"Cr"),
  std::pair<int, QString>(25,"Mn"),
  std::pair<int, QString>(26,"Fe"),
  std::pair<int, QString>(27,"Co"),
  std::pair<int, QString>(28,"Ni"),
  std::pair<int, QString>(29,"Cu"),
  std::pair<int, QString>(30,"Zn"),
  std::pair<int, QString>(31,"Ga"),
  std::pair<int, QString>(32,"Ge"),
  std::pair<int, QString>(33,"As"),
  std::pair<int, QString>(34,"Se"),
  std::pair<int, QString>(35,"Br"),
  std::pair<int, QString>(36,"Kr"),
  std::pair<int, QString>(37,"Rb"),
  std::pair<int, QString>(38,"Sr"),
  std::pair<int, QString>(39,"Y"),
  std::pair<int, QString>(40,"Zr"),
  std::pair<int, QString>(41,"Nb"),
  std::pair<int, QString>(42,"Mo"),
  std::pair<int, QString>(43,"Tc"),
  std::pair<int, QString>(44,"Ru"),
  std::pair<int, QString>(45,"Rh"),
  std::pair<int, QString>(46,"Pd"),
  std::pair<int, QString>(47,"Ag"),
  std::pair<int, QString>(48,"Cd"),
  std::pair<int, QString>(49,"In"),
  std::pair<int, QString>(50,"Sn"),
  std::pair<int, QString>(51,"Sb"),
  std::pair<int, QString>(52,"Te"),
  std::pair<int, QString>(53,"I"),
  std::pair<int, QString>(54,"Xe"),
  std::pair<int, QString>(55,"Cs"),
  std::pair<int, QString>(56,"Ba"),
  std::pair<int, QString>(57,"La"),
  std::pair<int, QString>(58,"Ce"),
  std::pair<int, QString>(59,"Pr"),
  std::pair<int, QString>(60,"Nd"),
  std::pair<int, QString>(61,"Pm"),
  std::pair<int, QString>(62,"Sm"),
  std::pair<int, QString>(63,"Eu"),
  std::pair<int, QString>(64,"Gd"),
  std::pair<int, QString>(65,"Tb"),
  std::pair<int, QString>(66,"Dy"),
  std::pair<int, QString>(67,"Ho"),
  std::pair<int, QString>(68,"Er"),
  std::pair<int, QString>(69,"Tm"),
  std::pair<int, QString>(70,"Yb"),
  std::pair<int, QString>(71,"Lu"),
  std::pair<int, QString>(72,"Hf"),
  std::pair<int, QString>(73,"Ta"),
  std::pair<int, QString>(74,"W"),
  std::pair<int, QString>(75,"Re"),
  std::pair<int, QString>(76,"Os"),
  std::pair<int, QString>(77,"Ir"),
  std::pair<int, QString>(78,"Pt"),
  std::pair<int, QString>(79,"Au"),
  std::pair<int, QString>(80,"Hg"),
  std::pair<int, QString>(81,"Tl"),
  std::pair<int, QString>(82,"Pb"),
  std::pair<int, QString>(83,"Bi"),
  std::pair<int, QString>(84,"Po"),
  std::pair<int, QString>(85,"At"),
  std::pair<int, QString>(86,"Rn"),
  std::pair<int, QString>(87,"Fr"),
  std::pair<int, QString>(88,"Ra"),
  std::pair<int, QString>(89,"Ac"),
  std::pair<int, QString>(90,"Th"),
  std::pair<int, QString>(91,"Pa"),
  std::pair<int, QString>(92,"U"),
  std::pair<int, QString>(93,"Np"),
  std::pair<int, QString>(94,"Pu"),
  std::pair<int, QString>(95,"Am"),
  std::pair<int, QString>(96,"Cm"),
  std::pair<int, QString>(97,"Bk"),
  std::pair<int, QString>(98,"Cf"),
  std::pair<int, QString>(99,"Es"),
  std::pair<int, QString>(100,"Fm"),
  std::pair<int, QString>(101,"Md"),
  std::pair<int, QString>(102,"No"),
  std::pair<int, QString>(103,"Lr"),
  std::pair<int, QString>(104,"Rf"),
  std::pair<int, QString>(105,"Db"),
  std::pair<int, QString>(106,"Sg"),
  std::pair<int, QString>(107,"Bh"),
  std::pair<int, QString>(108,"Hs"),
  std::pair<int, QString>(109,"Mt"),
  std::pair<int, QString>(110,"Ds"),
  std::pair<int, QString>(111,"Rg"),
  std::pair<int, QString>(112,"Cn")
};

static const std::map<int, QString> elementMap(elements,elements+sizeof(elements)/sizeof(std::pair<int,QString>));

#endif

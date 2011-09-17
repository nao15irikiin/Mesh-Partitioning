(* Recursive coordinate bisection is a technique for partitioning a list of coordinates
* into roughly equal sized partitions of spatially close elements. The general principle
* is find the dimension that spans the greatest distance and find a pivot that will
* split the list into two roughly equal halves. Then recurse down.
*)

(* Sort a list of floating point values *)
let sort_float_list l =
  List.sort (fun x y -> int_of_float (x -. y)) l

(* Find the span (max-min) of a list of floating point values *)
let span xs =
  let xs' = sort_float_list xs in
  let min = List.hd xs' in
  let max = List.nth xs' ((List.length xs') - 1) in
  max -. min

(* Find the midpoint of a list of floating point values
* nb: not so good since it doesn't take into account repeated values so the pivot
* won't necessarily divide the list into two equal halves.
*)
let get_pivot xs =
  let xs' = sort_float_list xs in
  let len = List.length xs' in
  if (len mod 2) = 0 then (*even*)
    let right = len / 2 in
    let left = right - 1 in
    ((List.nth xs' left) +. (List.nth xs' right)) /. 2.0
  else (*odd*)
    List.nth xs' (len/2)

exception NoLargestDim

(* Partition a list of (x,y,z) coordinates:
* - Find the dimension (x,y,z) with the largest span
* - Split this dimension into half by finding a pivot point
* - Recurse down until we get to the required depth
*
* We tag each coordinate (x,y,z) with a fourth component [n], the index value
* to form the 4-tuple (x,y,z,n). This is so we can extract the index value after
* partitioning to reorder other per-particle lists in the same order.
*)

let recursive_coordinate_bisection position depth =
  (* Extract the x,y,z component of a 4-tuple *)
  let get_x (x,_,_,_) = x in
  let get_y (_,y,_,_) = y in
  let get_z (_,_,z,_) = z in

  (* Tag each coordinate with its index value *)
  let append_index l =
    let rec f l n =
      match l with
      | [] -> []
      | (x,y,z)::xs -> (x,y,z,n) :: (f xs (n+1))
    in f l 0
  in
  let position' = append_index position in

  (* Recursive partition *)
  let rec partition_step position' depth =
    if depth = 0 then
      [position']
    else
      let xs = List.map get_x position' in
      let ys = List.map get_y position' in
      let zs = List.map get_z position' in
      let xs_span = span xs in
      let ys_span = span ys in
      let zs_span = span zs in
      (* Get pivot based on largest dimension *)
      let pivot, get_xyz, dim =
        if (xs_span >= ys_span && xs_span >= zs_span) then
          (get_pivot xs), get_x, "x"
        else if (ys_span >= xs_span && ys_span >= zs_span) then
          (get_pivot ys), get_y, "y"
        else if (zs_span >= xs_span && zs_span >= ys_span) then
          (get_pivot zs), get_z, "z"
        else
          let _ = Printf.printf "spans = {%f %f %f}\n" xs_span ys_span zs_span in
          raise NoLargestDim
      in
      (* Split the list according to the pivot *)
      let p1, p2 = List.partition (fun p -> (get_xyz p) <= pivot) position' in
      let _ =
        Printf.printf "Splitting %s dimension: %d/%d (%f)\n"
          dim (List.length p1) (List.length p2) pivot
      in
      (partition_step p1 (depth-1)) @ (partition_step p2 (depth-1))
  in

  partition_step position' depth



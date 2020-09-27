#include <igraph/igraph.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>

#define min(a, b) (((a) + (b) - fabs((a) - (b))) * 0.5)


void print_vector(igraph_vector_t *v) {
    long int i, l = igraph_vector_size(v);
    for (i = 0; i < l; i++) {
        printf(" %li", (long int) VECTOR(*v)[i]);
    }
    printf("\n");
}


void print_matrix(igraph_matrix_t *m) {
    long int i, j;
    for (i = 0; i < igraph_matrix_nrow(m); i++) {
        for (j = 0; j < igraph_matrix_ncol(m); j++) {
            printf(" %g", MATRIX(*m, i, j));
        }
        printf("\n");
    }
}


void print_vectorList(igraph_vector_ptr_t * vl) {
      long int i,n;
      n = igraph_vector_ptr_size(vl);
      for (i=0; i<n; i++) {
        igraph_vector_t * vt;
        vt=(igraph_vector_t*)igraph_vector_ptr_e(vl,i);
        print_vector((igraph_vector_t*)vt);
        igraph_vector_destroy(vt);
        free(vt);
      }
}


igraph_t build_auxiliary_node_connectivity(igraph_t *G)
{
	igraph_t H;
	igraph_vector_t edges, edges2, edges3, edges4;
	igraph_integer_t from, to, from2, to2;
	igraph_strvector_t id, id2;
	igraph_bool_t res;
	long int i, j; 
	long int nrow;
	long int ncol = 2;
	int G_vcount;
	char str[20];	char str2[20];
	char str3[20];	char str4[20];
	char A_name[20] = "A";
	char B_name[20] = "B";


	igraph_i_set_attribute_table(&igraph_cattribute_table);

	igraph_empty(&H, 0, 1);

	G_vcount = igraph_vcount(G);

	igraph_add_vertices(&H, G_vcount, 0);
	igraph_add_vertices(&H, G_vcount, 0);

	printf("H counts: %d\n", igraph_vcount(&H));


  	igraph_vector_init(&edges, 0);
  	for (i=0; i < igraph_vcount(G); i++) {
  		igraph_vector_push_back(&edges, i);
  		igraph_vector_push_back(&edges, igraph_vcount(G) + i);
  	}
  
   	igraph_add_edges(&H, &edges, 0);
   	igraph_vector_destroy(&edges);


   	igraph_vector_init(&edges2, 0);
   	igraph_vector_init(&edges3, 0);
   	for (i = 0; i < igraph_ecount(G); i++) {
   		igraph_edge(G, i, &from, &to);

   		igraph_vector_push_back(&edges2, to+igraph_vcount(G));
   		igraph_vector_push_back(&edges2, from);
   		igraph_vector_push_back(&edges3, from+igraph_vcount(G));
   		igraph_vector_push_back(&edges3, to);
   	}

   	igraph_add_edges(&H, &edges2, 0);
  	igraph_vector_destroy(&edges2);

    igraph_add_edges(&H, &edges3, 0);
  	igraph_vector_destroy(&edges3);



  	igraph_vector_init(&edges4, igraph_ecount(&H));
  	igraph_vector_fill(&edges4, 1);
 	SETEANV(&H, "capacity", &edges4);
 	igraph_vector_destroy(&edges4);

    return H;
}


igraph_t build_residual_network(igraph_t *H)
{
	igraph_integer_t from, to, from_2, to_2, from_3, to_3;
	igraph_vector_t edges, edges2, edges3;
	igraph_t R;
	int kbs;
	long int i;


	igraph_empty(&R, 0, 1);

	igraph_add_vertices(&R, igraph_vcount(H), 0);

   	igraph_vector_init(&edges, 0);
   	igraph_vector_init(&edges2, 0);
   	for (i = 0; i < igraph_ecount(H); i++) {
   		igraph_edge(H, i, &from, &to);
   		igraph_vector_push_back(&edges, from);
   		igraph_vector_push_back(&edges, to);
   		igraph_vector_push_back(&edges2, to);
   		igraph_vector_push_back(&edges2, from);
   	}

  	igraph_add_edges(&R, &edges, 0);
  	igraph_vector_destroy(&edges);

 	igraph_vector_init(&edges, (int)igraph_ecount(&R));
  	igraph_vector_fill(&edges, 1);
  	SETEANV(&R, "capacity", &edges);
 	igraph_vector_destroy(&edges);

 	igraph_add_edges(&R, &edges2, 0);
 	igraph_vector_init(&edges3, (int)igraph_vector_size(&edges2));
  	igraph_vector_fill(&edges3, 0);

  	for (i = (int)igraph_ecount(&R)-(int)igraph_vector_size(&edges3)/2; i < (int)igraph_ecount(&R); i++) {
  		SETEAN(&R, "capacity", i, 0);
  	}
  	igraph_vector_destroy(&edges3);

 	return R;
}


void k_max_ids(igraph_t *G, long int k, igraph_vector_t *k_vec_r, igraph_vector_t *k_vec_ids_r)
{
	igraph_vector_t result, edges, k_vec, k_vec_ids;
	igraph_matrix_t m;
	igraph_vit_t vit;
	igraph_vs_t vs;
	long int i, j;


	igraph_vector_init(&result, igraph_vcount(G));
	igraph_degree(G, &result, igraph_vss_all(), IGRAPH_ALL, IGRAPH_NO_LOOPS);

	igraph_matrix_init(&m, igraph_vcount(G), 2);

	igraph_vs_seq(&vs, 0, igraph_vcount(G));
	igraph_vit_create(G, vs, &vit);

	igraph_vector_init(&edges, igraph_vcount(G));
	for (i = 0; i < igraph_vcount(G); i++) {
		VECTOR(edges)[i] = (long int)IGRAPH_VIT_GET(vit);
		IGRAPH_VIT_NEXT(vit);
	}

	igraph_matrix_set_col(&m, &edges, 0);
	igraph_matrix_set_col(&m, &result, 1);

	igraph_vector_sort(&result);

	igraph_vector_init(&k_vec, k);
	igraph_vector_init(&k_vec_ids, k);
	for (j = 1; j < k+1; j++) {
		for (i = 0; i < igraph_vcount(G); i++) {
			if ((int)igraph_matrix_e(&m, i, 1) == (int)VECTOR(result)[igraph_vcount(G)-j]) {
				if (j <= 1) {
					VECTOR(k_vec)[j-1] = (int)igraph_matrix_e(&m, i, 1);
					VECTOR(k_vec_ids)[j-1] = (int)igraph_matrix_e(&m, i, 0);
				} else {
					if (VECTOR(k_vec_ids)[j-2] != (int)igraph_matrix_e(&m, i, 0)) {
						VECTOR(k_vec)[j-1] = (int)igraph_matrix_e(&m, i, 1);
						VECTOR(k_vec_ids)[j-1] = (int)igraph_matrix_e(&m, i, 0);						
					}
				}				
			}
		}
	}
	*k_vec_r = k_vec;
	*k_vec_ids_r = k_vec_ids;

	return;
}


int is_separatings_set(igraph_t *G, igraph_vector_t *k_vec_ids_r)
{
	igraph_t S;
	igraph_vs_t vs;
	igraph_bool_t res;


	if (igraph_vector_size(k_vec_ids_r) == igraph_vcount(G) - 1) {
		return 1;
	}

	igraph_copy(&S, G);
	igraph_vs_vector(&vs, k_vec_ids_r);
	igraph_delete_vertices(&S, vs);
	igraph_is_connected(&S, &res, IGRAPH_WEAK);
	
	if (res == 0) {
		return 1;
	} else if (res == 1) {
		return 0;
	}

}


igraph_vector_t distinct_elements(igraph_vector_t *vt) {

	int i, j, n;
	igraph_vector_t res;

	n = igraph_vector_size(vt);

	igraph_vector_init(&res, 0);
	for (i=0; i<n; i++) {
		for (j=0; j<i; j++) {
			if ((int)VECTOR(*vt)[i] == (int)VECTOR(*vt)[j])
				break;
		}
		if (i == j)
			igraph_vector_push_back(&res, (int)VECTOR(*vt)[i]);
	}
	return res;
}


int find_elem(int elem, igraph_vector_t *array) {

	int i;
	int size;

	size = igraph_vector_size(array);

   	for (i = 0; i < size; i++) {

    	if ((int)VECTOR(*array)[i] == elem) {
     		return 1;
       	}
   	}

   	if ((int)VECTOR(*array)[i] != elem) {
       	return 0;
   	}
  	return 0;
}


igraph_es_t pair_u_v(int u, int v) {
	igraph_vector_t pair;
	igraph_es_t es;

	igraph_vector_init(&pair, 0);
	igraph_vector_push_back(&pair, u);
	igraph_vector_push_back(&pair, v);

	igraph_es_pairs(&es, &pair, 1);

	igraph_vector_destroy(&pair);

	return es;
}


int edmonds_karp_core(igraph_t *R, igraph_t *G, int s, int t, int cutoff, igraph_t *R_result) {
	int inf, flow_value, u_w;
	igraph_vector_t path;
	int res_v;
	igraph_vector_t res_pred_dict_all_keys;
	igraph_vector_t res_succ_dict_all_keys;
	igraph_vector_t res_pred_dict_all_values;
	igraph_vector_t res_succ_dict_all_values;
	long int i;


	inf = 144;

	int augment(igraph_vector_t *path, igraph_t *R_a, igraph_t *R_res_a) {
		igraph_eit_t it_3, it_4, it_5;
		igraph_es_t es_3, es_4, es_5;
		int u, flow, flow_capa_id, attr_capa, attr_flow, pair_id_1, pair_id_2;
		long int v;

		R_a = R;

		flow = inf;

		u = (int)VECTOR(*path)[0];

		for (v = 1; v < igraph_vector_size(path); v++) {

			es_3 = pair_u_v(u, (int)VECTOR(*path)[v]);

			igraph_eit_create(R_a, es_3, &it_3);

			flow_capa_id = (int)IGRAPH_EIT_GET(it_3);
			
			attr_capa = (int)EAN(R_a, "capacity", flow_capa_id);
			attr_flow = (int)EAN(R_a, "flow", flow_capa_id);
			flow = min(flow, attr_capa-attr_flow);

			u = (int)VECTOR(*path)[v];
		}

		u = (int)VECTOR(*path)[0];

		for (v = 1; v < igraph_vector_size(path); v++) {

			es_4 = pair_u_v(u, (int)VECTOR(*path)[v]);

			igraph_eit_create(R_a, es_4, &it_4);

			pair_id_1 = (int)IGRAPH_EIT_GET(it_4);

			SETEAN(R, "flow", pair_id_1, EAN(R, "flow", pair_id_1) + flow);

			es_5 = pair_u_v(u, (int)VECTOR(*path)[v]);

			igraph_eit_create(R_a, es_5, &it_5);

			pair_id_2 = (int)IGRAPH_EIT_GET(it_5);

			SETEAN(R_a, "flow", pair_id_2, EAN(R_a, "flow", pair_id_2) + flow);

			u = (int)VECTOR(*path)[v];

		}

		*R_res_a = *R_a;

		return flow;
	}


	int bidirectional_bfs(igraph_t *G_b, igraph_t *R_b, int *res_v, igraph_vector_t *res_pred_dict_all_keys, igraph_vector_t *res_succ_dict_all_keys, igraph_vector_t *res_pred_dict_all_values, igraph_vector_t *res_succ_dict_all_values) {
		igraph_vector_t q_s, pred_dict_all_keys, q_t, succ_dict_all_keys;
		igraph_vector_t q_vec, restricted_1, order_1, rank_1, father_1, pred_1, succ_1, dist_1, restricted_2, order_2, rank_2, father_2, pred_2, succ_2, dist_2;
		igraph_vector_t R_succ, R_pred, q_ss, q_tt;
		igraph_es_t es_1, es_2;
		igraph_eit_t it_1, it_2;

		char dict_key_1, dict_key_2, R_succ_v_char, R_pred_v_char;
		int dict_value_1, dict_value_2, R_succ_v, R_pred_v, elem_1, elem_2, elem_3, elem_4, capa_val, flow_val, flow_value, flow_capa_id, capa_val_2, flow_val_2, flow_capa_id_2;
		long int u, v, w, x, o, p, i;

		igraph_vector_init(&q_s, 1);
		igraph_vector_fill(&q_s, s + igraph_vcount(G_b));
		igraph_vector_init(&pred_dict_all_keys, 0);
		
		dict_key_1 = s + igraph_vcount(G_b);
		dict_value_1 = -9;

		igraph_vector_t pred_dict_all_values;

		igraph_vector_init(&pred_dict_all_keys, 0);
		igraph_vector_init(&pred_dict_all_values, 0);

		igraph_vector_push_back(&pred_dict_all_keys, dict_key_1);
		igraph_vector_push_back(&pred_dict_all_values, dict_value_1);

		igraph_vector_init(&q_t, 1);
		igraph_vector_fill(&q_t, t);
		igraph_vector_init(&succ_dict_all_keys, 0);

		dict_key_2 = t;
		dict_value_2 = -9;

		igraph_vector_t succ_dict_all_values;

		igraph_vector_init(&succ_dict_all_keys, 0);
		igraph_vector_init(&succ_dict_all_values, 0);

		igraph_vector_push_back(&succ_dict_all_keys, dict_key_2);
		igraph_vector_push_back(&succ_dict_all_values, dict_value_2);

		while (1) {
			igraph_vector_init(&q_vec, 0);
			if (igraph_vector_size(&q_s) <= igraph_vector_size(&q_t)) {
				
				for (u = 0; u < igraph_vector_size(&q_s); u++) {

					igraph_vector_init(&restricted_1, 0);
					igraph_vector_init(&order_1, 0);
					igraph_vector_init(&rank_1, 0);
					igraph_vector_init(&father_1, 0);
					igraph_vector_init(&pred_1, 0);
					igraph_vector_init(&succ_1, 0);
					igraph_vector_init(&dist_1, 0);

					for (o = 0; o < igraph_vcount(R_b); o++) {
						igraph_vector_push_back(&restricted_1, o);
					}

					igraph_vector_init(&q_ss, 1);
					VECTOR(q_ss)[0] = (int)VECTOR(q_s)[u];

					igraph_bfs(R_b, (int)VECTOR(q_s)[u], &q_ss, IGRAPH_ALL, 0, &restricted_1,
						&order_1, &rank_1, &father_1, &pred_1, &succ_1, &dist_1, 0, 0);

					igraph_vector_init(&R_succ, 0);
					for (p = 0; p < igraph_vector_size(&dist_1); p++) {
						if ((int)VECTOR(dist_1)[p] == 1) {
							igraph_vector_push_back(&R_succ, p);
						}
					}

					igraph_vector_destroy(&restricted_1);
					igraph_vector_destroy(&order_1);
					igraph_vector_destroy(&rank_1);
					igraph_vector_destroy(&father_1);
					igraph_vector_destroy(&pred_1);
					igraph_vector_destroy(&succ_1);
					igraph_vector_destroy(&dist_1);

					for (v = 0; v < igraph_vector_size(&R_succ); v++) {

						R_succ_v = (int)VECTOR(R_succ)[v];

						elem_1 = find_elem(R_succ_v, &pred_dict_all_keys);
						elem_2 = find_elem(R_succ_v, &succ_dict_all_keys);

						es_1 = pair_u_v((int)VECTOR(q_s)[u], R_succ_v);

						igraph_eit_create(R_b, es_1, &it_1);

						flow_capa_id = (int)IGRAPH_EIT_GET(it_1);

						//igraph_es_destroy(&es_1);
						//igraph_eit_destroy(&it_1);
						
						capa_val = (int)EAN(R_b, "capacity", flow_capa_id);
						flow_val = (int)EAN(R_b, "flow", flow_capa_id);

						if (elem_1 != 1 && flow_val < capa_val) {
							igraph_vector_push_back(&pred_dict_all_keys, R_succ_v);
							igraph_vector_push_back(&pred_dict_all_values, (int)VECTOR(q_s)[u]);
							if (elem_2 == 1) {
								*res_v = R_succ_v;
								*res_pred_dict_all_keys = pred_dict_all_keys;
								*res_succ_dict_all_keys = succ_dict_all_keys;
								*res_pred_dict_all_values = pred_dict_all_values;
								*res_succ_dict_all_values = succ_dict_all_values;
								return 0;
							}

							igraph_vector_push_back(&q_vec, R_succ_v);
						}
					}

				}
				if ((int)igraph_vector_empty(&q_vec) == 1) {
					*res_v = -9;
					return 0;
				}
				igraph_vector_destroy(&q_s);
				igraph_vector_copy(&q_s, &q_vec);

			} else {

				for (w = 0; w < igraph_vector_size(&q_t); w++) {

					igraph_vector_init(&restricted_2, 0);
					igraph_vector_init(&order_2, 0);
					igraph_vector_init(&rank_2, 0);
					igraph_vector_init(&father_2, 0);
					igraph_vector_init(&pred_2, 0);
					igraph_vector_init(&succ_2, 0);
					igraph_vector_init(&dist_2, 0);

					for (o = 0; o < igraph_vcount(R_b); o++) {
						igraph_vector_push_back(&restricted_2, o);
					}

					igraph_vector_init(&q_tt, 1);
					VECTOR(q_tt)[0] = (int)VECTOR(q_t)[w];

					igraph_bfs(R_b, (int)VECTOR(q_t)[w], &q_tt, IGRAPH_ALL, 0, &restricted_2,
						&order_2, &rank_2, &father_2, &pred_2, &succ_2, &dist_2, 0, 0);

					igraph_vector_init(&R_pred, 0);
					for (p = 0; p < igraph_vector_size(&dist_2); p++) {
						if ((int)VECTOR(dist_2)[p] == 1) {
							igraph_vector_push_back(&R_pred, p);
						}
					}
					igraph_vector_destroy(&restricted_2);
					igraph_vector_destroy(&order_2);
					igraph_vector_destroy(&rank_2);
					igraph_vector_destroy(&father_2);
					igraph_vector_destroy(&pred_2);
					igraph_vector_destroy(&succ_2);
					igraph_vector_destroy(&dist_2);
				
					for (x = 0; x < igraph_vector_size(&R_pred); x++) {
						R_pred_v = (int)VECTOR(R_pred)[x];

						elem_3 = find_elem(R_pred_v, &succ_dict_all_keys);
						elem_4 = find_elem(R_pred_v, &pred_dict_all_keys);
						
						es_2 = pair_u_v(R_pred_v, (int)VECTOR(q_t)[w]);

						igraph_eit_create(R_b, es_2, &it_2);

						flow_capa_id_2 = (int)IGRAPH_EIT_GET(it_2);

						//igraph_es_destroy(&es_2);
						//igraph_eit_destroy(&it_2);
						
						capa_val_2 = (int)EAN(R_b, "capacity", flow_capa_id_2);
						flow_val_2 = (int)EAN(R_b, "flow", flow_capa_id_2);
						// printf("\n\n");
						if (elem_3 != 1 && flow_val_2 < capa_val_2) {
							igraph_vector_push_back(&succ_dict_all_keys, R_pred_v);
							igraph_vector_push_back(&succ_dict_all_values, (int)VECTOR(q_t)[w]);
							if (elem_4 == 1) {
								*res_v = R_pred_v;
								*res_pred_dict_all_keys = pred_dict_all_keys;
								*res_succ_dict_all_keys = succ_dict_all_keys;
								*res_pred_dict_all_values = pred_dict_all_values;
								*res_succ_dict_all_values = succ_dict_all_values;
								return 0;
							}
							igraph_vector_push_back(&q_vec, R_pred_v);
						}
					}
				}
				if ((int)igraph_vector_empty(&q_vec) == 1) {
					*res_v = -9;
					return 0;
				}
				igraph_vector_destroy(&q_t);
				igraph_vector_copy(&q_t, &q_vec);

			}

			igraph_vector_destroy(&q_vec);
		}
	}


	flow_value = 0;

	while (flow_value < cutoff) {
		igraph_t R_res;
		igraph_vector_t path_w;
		int kbs;

		kbs = bidirectional_bfs(G, R, &res_v, &res_pred_dict_all_keys, &res_succ_dict_all_keys, &res_pred_dict_all_values, &res_succ_dict_all_values);

		if (res_v == -9) {
			printf("ddd\n");
			break;
		}

		igraph_vector_init(&path_w, 0);
		igraph_vector_push_back(&path_w, res_v);

		u_w = res_v;

		while (u_w != s + igraph_vcount(G)) {
			for (i = 0; i < igraph_vector_size(&res_pred_dict_all_keys); i++) {
				if ((int)VECTOR(res_pred_dict_all_keys)[i] == u_w) {
					u_w = (int)VECTOR(res_pred_dict_all_values)[i];
				}
			}
			igraph_vector_push_back(&path_w, u_w);
			print_vector(&path_w);
		}

		igraph_vector_reverse(&path_w);

		print_vector(&path_w);

		u_w = res_v;

		while (u_w != t) {
			for (i = 0; i < igraph_vector_size(&res_succ_dict_all_keys); i++) {
				if ((int)VECTOR(res_succ_dict_all_keys)[i] == u_w) {
					u_w = (int)VECTOR(res_succ_dict_all_values)[i];
				}
			}
			igraph_vector_push_back(&path_w, u_w);
		}

		flow_value += augment(&path_w, R, &R_res);

		R = &R_res;

	}

	igraph_vector_destroy(&res_pred_dict_all_keys);
	igraph_vector_destroy(&res_pred_dict_all_values);
	igraph_vector_destroy(&res_succ_dict_all_keys);
	igraph_vector_destroy(&res_succ_dict_all_values);


	*R_result = *R;

	return flow_value;
}


igraph_t condensation(igraph_t *R_res, igraph_vector_t *c_mapping_key, igraph_vector_t *c_mapping_value) {
	igraph_vector_t membership, csize, mapping_key, mapping_value, mem_edges;
	igraph_t C;
	int no, from, to;
	long int i;

	igraph_vector_init(&membership, 0);
   	igraph_vector_init(&csize, 0);

   	igraph_clusters(R_res, &membership, &csize, &no, IGRAPH_STRONG);

   	print_vector(&membership);
   	print_vector(&csize);
   	printf("%d\n", (int)no);

	igraph_empty(&C, 0, 1);

	igraph_add_vertices(&C, (int)no, 0);

	igraph_vector_init(&mapping_key, 0);
	igraph_vector_init(&mapping_value, 0);

	for (i = 0; i < igraph_vector_size(&membership); i++) {
		igraph_vector_push_back(&mapping_key, i);
		igraph_vector_push_back(&mapping_value, VECTOR(membership)[i]);
	}

	igraph_vector_init(&mem_edges, 0);

	for (i = 0; i < igraph_ecount(R_res); i++) {
		igraph_edge(R_res, i, &from, &no);

		if (VECTOR(mapping_value)[from] != VECTOR(mapping_value)[no]) {
			igraph_vector_push_back(&mem_edges, VECTOR(mapping_value)[from]);
			igraph_vector_push_back(&mem_edges, VECTOR(mapping_value)[no]);
		}
	}

	igraph_add_edges(&C, &mem_edges, 0);

	*c_mapping_key = mapping_key;
	*c_mapping_value = mapping_value;

	return C;
}


igraph_vector_t descendants_at_distance(igraph_t *G, int source, int distance) {
	int current_distance, vertex, child;
	igraph_vector_t queue, visited, next_vertices;
	igraph_vector_t set, eids;

	igraph_vector_init(&set, 0);

	current_distance = 0;

	igraph_vector_init(&queue, 0);
	igraph_vector_push_back(&queue, source);

	igraph_vector_init(&visited, 0);
	igraph_vector_push_back(&visited, source);

	while (igraph_vector_empty(&queue) == 0) {

		if (current_distance == distance) {
			return queue;
		}

		current_distance += 1;

		igraph_vector_init(&next_vertices, 0);

		for (vertex = 0; vertex < igraph_vector_size(&queue); vertex++) {
			igraph_vector_init(&eids, 0);
			igraph_incident(G, &eids, VECTOR(queue)[vertex], IGRAPH_OUT);

			if (igraph_vector_empty(&eids) == 0) {
				for (child = 0; child < igraph_vector_size(&eids); child++) {
					if (find_elem(VECTOR(eids)[child], &visited) == 0) {
						igraph_vector_push_back(&visited, VECTOR(eids)[child]);
						igraph_vector_push_back(&next_vertices, VECTOR(eids)[child]);
					}
				}
			}	
		}

		igraph_vector_copy(&queue, &next_vertices);
	}
	return set;
}


igraph_t transitive_closure_dag(igraph_t *L, igraph_vector_t *topo_order) {
	igraph_t TC;
	igraph_vector_t dad, TC_edges;
	long int v, u;

	igraph_copy(&TC, L);

	igraph_vector_reverse(topo_order);

	for (v = 0; v < igraph_vector_size(topo_order); v++) {
		dad = descendants_at_distance(&TC, v, 2);

		if (igraph_vector_empty(&dad) == 0) {
			for (u = 0; u < igraph_vector_size(&dad); u++) {
				igraph_vector_init(&TC_edges, 0);
				igraph_vector_push_back(&TC_edges, u);
				igraph_vector_push_back(&TC_edges, v);
				igraph_vector_destroy(&TC_edges);
			}
			igraph_add_edges(&TC, &TC_edges, 0);
		}

	}
	return TC;
}


igraph_vector_ptr_t antichains(igraph_t *L) {
	igraph_vector_ptr_t antichains, stacks, resulting_antichains;
	igraph_vector_t *current_antichain, *current_stack, *new_antichain, *new_stack, *new_antichain_w, *new_stack_w, tcx_neis, tct_neis;
	igraph_integer_t x_w;
	igraph_t TC;
	long int t;

	igraph_vector_ptr_init(&antichains, 0);
	igraph_vector_ptr_init(&stacks, 0);
	igraph_vector_ptr_init(&resulting_antichains, 0);

	new_stack = malloc(sizeof(igraph_vector_t));
	igraph_vector_init(new_stack, 0);
	igraph_topological_sorting(L, new_stack, IGRAPH_OUT);
	TC = transitive_closure_dag(L, new_stack);

	igraph_vector_ptr_push_back(&stacks, new_stack);

	new_antichain = malloc(sizeof(igraph_vector_t));
	igraph_vector_init(new_antichain, 0);
	igraph_vector_ptr_push_back(&antichains, new_antichain);
	//igraph_vector_ptr_push_back(&antichains, NULL);

	while (igraph_vector_ptr_empty(&stacks) == 0) {

		//if (igraph_vector_ptr_size(&antichains) > 0) {
		current_antichain = igraph_vector_ptr_pop_back(&antichains);
		igraph_vector_ptr_push_back(&resulting_antichains, current_antichain);
		//}
		current_stack = igraph_vector_ptr_pop_back(&stacks);

		printf("current_stack:\n");
		print_vector(current_stack);
		printf("current_antichain:\n");
		print_vector(current_antichain);

		while (igraph_vector_empty(current_stack) == 0) {
			new_antichain_w = malloc(sizeof(igraph_vector_t));
			new_stack_w = malloc(sizeof(igraph_vector_t));

			igraph_vector_init(new_antichain_w, 0);
			igraph_vector_init(new_stack_w, 0);

			//if (igraph_vector_ptr_size(&antichains) > 0) {
			igraph_vector_copy(new_antichain_w, current_antichain);
			//}

			x_w = igraph_vector_pop_back(current_stack);
			igraph_vector_push_back(new_antichain_w, x_w);

   			if (igraph_vector_empty(current_stack) == 0) {
   				for (t = 0; t < igraph_vector_size(current_stack); t++) {
   					igraph_vector_init(&tcx_neis, 0);
   					igraph_vector_init(&tct_neis, 0);

  					igraph_neighbors(&TC, &tcx_neis, x_w, IGRAPH_OUT);
  					igraph_neighbors(&TC, &tct_neis, (int)VECTOR(*current_stack)[t], IGRAPH_OUT);

					if (find_elem(VECTOR(*current_stack)[t], &tcx_neis) == 0 && find_elem(x_w, &tct_neis) == 0) {
						igraph_vector_push_back(new_stack, VECTOR(*current_stack)[t]);
					}
   				}
   			}
   			printf("new_antichain_w:\n");
   			print_vector(new_antichain_w);
   			printf("new_stack_w:\n");
   			print_vector(new_stack_w);

			igraph_vector_ptr_push_back(&antichains, new_antichain_w);
			igraph_vector_ptr_push_back(&stacks, new_stack_w);

			igraph_vector_copy(current_antichain, new_antichain_w);

		   	printf("current_antichain:\n");
   			print_vector(current_antichain);

		}

		igraph_vector_destroy(current_stack);
		igraph_free(current_stack);
	}

	igraph_vector_ptr_destroy(&stacks);
	igraph_vector_ptr_destroy(&antichains);


	printf("size: %d\n", igraph_vector_ptr_size(&resulting_antichains));

	return resulting_antichains;
}


int main(void)
{
	igraph_t G;
	igraph_t H;
	igraph_t R;
	igraph_t G_r;
	igraph_t R_res;
	igraph_t L;
	FILE *file;
	long int k, sep;
	igraph_vector_t k_vec_r, k_vec_ids_r, ids;

	igraph_vector_ptr_t res;
	igraph_vs_t vids;
	long int i, j, l, m, n, q, G_vcount;
	igraph_vector_t *vt;
	igraph_vector_t uniq_res, edges;
	igraph_vs_t uniq_vs;
	igraph_matrix_t m_r;

	igraph_integer_t size;
	long int attr_capa, attr_flow;
	int pair_id_1, pair_id_2, cutoff;

	char mat_e[20];
	char data[G_vcount][2][20];
	int u_w, j_node, k_node;

	igraph_vector_t q_s, pred_dict_all_keys, q_t, succ_dict_all_keys;
	igraph_vector_t q_vec, restricted_1, order_1, rank_1, father_1, pred_1, succ_1, dist_1, restricted_2, order_2, rank_2, father_2, pred_2, succ_2, dist_2;
	igraph_vector_t R_succ, R_pred, empty_vector, empty_vector_2, empty_vector_3, empty_vector_4;
	igraph_es_t es_1, es_2, es_3, es_4, es_5;
	igraph_eit_t it_1, it_2;

	int dict_key_1, dict_value_1, dict_key_2, dict_value_2, R_succ_v, R_pred_v, elem_1, elem_2, capa_val, flow_val, flow_value, flow_capa_id, empty_v, empty_v_2;
	int s, t, nn, scc;
	long int u, v, w, x, o, p, k2, b, c;

	igraph_vector_t saturated_edges, saturated_edges_ids, S;
	int from, to, from_2, to_2, capa_val_sat, flow_val_sat;
	igraph_es_t es_sat, es_sat_del;
	igraph_eit_t it_sat;

	igraph_vector_ptr_t antichains_res, cutset;

	igraph_vector_t cmap_key, cmap_value, hu_neis, *u_w_pair, *temp_edge, node_cut, node_cut_uniq;


	igraph_i_set_attribute_table(&igraph_cattribute_table);

	file = fopen("F:/practice/k_components/cornwell_burchard_2019.txt", "r");

	if (!file) {
		return 1;
	}

	igraph_read_graph_edgelist(&G, file, 0, 1);

	fclose(file);


	G_vcount = igraph_vcount(&G);

	H = build_auxiliary_node_connectivity(&G);

	R = build_residual_network(&H);

	k = 3;

	k_max_ids(&G, k, &k_vec_r, &k_vec_ids_r);

	sep = is_separatings_set(&G, &k_vec_ids_r);


	printf("k_vec_ids_r size: %d\n", igraph_vector_size(&k_vec_ids_r));
	print_vector(&k_vec_ids_r);

	for (x = 0; x < igraph_vector_size(&k_vec_ids_r); x++) {

		igraph_vs_vector(&vids, &k_vec_ids_r);
		igraph_vector_ptr_init(&res, 0);
		igraph_neighborhood(&G, &res, vids, 1, IGRAPH_ALL, 1);

		vt = (igraph_vector_t*)igraph_vector_ptr_e(&res, x);

		for (i = 0; i < igraph_vector_size(&k_vec_ids_r); i++) {
			igraph_vector_push_back(vt, VECTOR(k_vec_ids_r)[i]);
		}
		
		uniq_res = distinct_elements(vt);

		igraph_copy(&G_r, &G);

		for (j = 0; j < igraph_vcount(&G_r); j++) {
			SETVAN(&G_r, "prev_ids", j, j);
		}

		igraph_vs_vector(&uniq_vs, &uniq_res);
		igraph_delete_vertices(&G_r, uniq_vs);

		//v = 1;
		for (v = 0; v < igraph_vcount(&G_r); v++) {
			j_node = (int)VECTOR(k_vec_ids_r)[x];
			k_node = (int)VAN(&G_r, "prev_ids", v);

			igraph_vector_init(&edges, igraph_ecount(&R));
			igraph_vector_fill(&edges, 0);
			SETEANV(&R, "flow", &edges);

			cutoff = 100;

			flow_value = edmonds_karp_core(&R, &G, j_node, k_node, cutoff, &R_res);

			if (flow_value == k) {
				igraph_vector_init(&saturated_edges, 0);

			   	for (b = 0; b < igraph_ecount(&R_res); b++) {
			   		igraph_edge(&R_res, b, &from, &to);

			   		es_sat = pair_u_v(from, to);

			   		igraph_eit_create(&R_res, es_sat, &it_sat);

					flow_capa_id = (int)IGRAPH_EIT_GET(it_sat);

					//igraph_es_destroy(&es_sat);
					//igraph_eit_destroy(&it_sat);
			
					capa_val_sat = (int)EAN(&R_res, "capacity", flow_capa_id);
					flow_val_sat = (int)EAN(&R_res, "flow", flow_capa_id);

			   	}

			   	igraph_es_pairs(&es_sat_del, &saturated_edges, 1);

				igraph_delete_edges(&R_res, es_sat_del);


		 	   	L = condensation(&R_res, &cmap_key, &cmap_value);

		 	   	print_vector(&cmap_key);
		 	   	print_vector(&cmap_value);


			   	antichains_res = antichains(&L);

			   	igraph_vector_init(&S, 0);

			   	for (i = 0; i < igraph_vector_ptr_size(&antichains_res); i++) {
			   		for (j = 0; j < igraph_vector_size(&cmap_key); j++) {
			   			nn = VECTOR(cmap_key)[j];
			   			scc = VECTOR(cmap_value)[j];
			   			if (find_elem(scc, igraph_vector_ptr_e(&antichains_res, i)) == 1) {
			   				if (find_elem(nn, &S) == 0) {
			   					igraph_vector_push_back(&S, nn);
			   				}
			   			}
			   		}
		   		
			   		igraph_vector_ptr_init(&cutset, 0);
			   		for (i = 0; i < igraph_vector_size(&S); i++) {
			   		   	igraph_vector_init(&hu_neis, 0);
			   			igraph_neighbors(&H, &hu_neis, (int)VECTOR(S)[i], IGRAPH_OUT);
			   			for (w = 0; w < igraph_vector_size(&hu_neis); w++) {
			   				if (find_elem((int)VECTOR(hu_neis)[w], &S) == 0) {
			   					u_w_pair = malloc(sizeof(igraph_vector_t));
								igraph_vector_init(u_w_pair, 0);
								igraph_vector_push_back(u_w_pair, (int)VECTOR(S)[i]);
								igraph_vector_push_back(u_w_pair, (int)VECTOR(hu_neis)[w]);				
								igraph_vector_ptr_push_back(&cutset, u_w_pair);
			   				} else {
			   					u_w_pair = malloc(sizeof(igraph_vector_t));
								igraph_vector_init(u_w_pair, 0);
								igraph_vector_ptr_push_back(&cutset, u_w_pair);
			   				}
			   			}
			   			igraph_vector_destroy(&hu_neis);
			   		}


		   			igraph_vector_init(&node_cut, 0);
					for (i = 0; i < igraph_vector_ptr_size(&cutset); i++) {
			   			temp_edge = igraph_vector_ptr_e(&cutset, i);
			   			for (j = 0; j < igraph_vector_size(temp_edge); j++) {
			   				if(find_elem((int)VECTOR(*temp_edge)[j], &node_cut) == 0) {
			   					if ((int)VECTOR(*temp_edge)[j] < igraph_vcount(&G)) {
			   						igraph_vector_push_back(&node_cut, (int)VECTOR(*temp_edge)[j]);
			   					} else {
			   						igraph_vector_push_back(&node_cut, (int)VECTOR(*temp_edge)[j] - igraph_vcount(&G));
			   					}
			   				}
			   			}
			   		}
			   		node_cut_uniq = distinct_elements(&node_cut);

			   		print_vector(&node_cut_uniq);

			   		if (igraph_vector_size(&node_cut_uniq) == k) {
			   			printf("node_cut: \n");
			   			print_vector(&node_cut_uniq);
					}
				}
			}
		}
	}

	return 0;
}
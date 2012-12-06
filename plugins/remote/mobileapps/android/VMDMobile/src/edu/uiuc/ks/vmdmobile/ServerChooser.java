package edu.uiuc.ks.vmdmobile;
 
import android.app.ListActivity;
import android.content.Context;
import android.content.SharedPreferences;
import android.os.Bundle;
import android.preference.Preference;
import android.preference.PreferenceActivity;
import android.preference.Preference.OnPreferenceClickListener;
import android.view.View;
import android.widget.ArrayAdapter;
import android.widget.ListView;
import android.widget.Toast;
 
public class ServerChooser extends ListActivity {
   @Override
   protected void onCreate(Bundle savedInstanceState) {
      super.onCreate(savedInstanceState);
      String[] names = new String[] { "Add New VMD Server...",
         "athine.ks.uiuc.edu"};
      // Create an ArrayAdapter, that will actually make the Strings above
      // appear in the ListView
      this.setListAdapter(new ArrayAdapter<String>(this,
                                          android.R.layout.simple_list_item_1,
                                          names));
   }

   @Override
   protected void onListItemClick(ListView l, View v, int position, long id) {
      super.onListItemClick(l, v, position, id);
      // Get the item that was clicked
      Object o = this.getListAdapter().getItem(position);
      String keyword = o.toString();
      Toast.makeText(this, "You selected " + position + ": " + keyword, Toast.LENGTH_LONG)
            .show();
   }

}
